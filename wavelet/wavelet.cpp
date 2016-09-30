#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_io.hpp>
using namespace boost::gil;

#include "..\gil_utils\color_arithm.h"
#include "..\gil_utils\float_views_io.h"

// заполняет в каждой строке писелы [0,x0), копируя их с зеркальной симметрией из пикселов [x0,...)
template <typename V>
void fill_left_border(const V & view, int x0)
{
  if (2 * x0 > view.width())
    throw std::runtime_error("too narrow view given in fill_left_border");
  for (int y = 0; y < view.height(); ++y)
  {
    const auto row = view.row_begin(y);
    for (int x = 0; x < x0; ++x)
      *(row + x) = *(row + 2 * x0 - x - 1);
  }
}

// заполняет в каждой строке первые brd-пикселов и последние brd-пикселов, копируя их с зеркальной симметрией их пикселов плиже к центру;
// то же делает и в каждом столбце;
// значения пикселов в углах изображения размером brd x brd меняется на неопределенное
template <typename V>
void fill_all_borders(const V & view, int brd)
{
  fill_left_border(view, brd);
  fill_left_border(rotated180_view(view), brd);
  fill_left_border(rotated90cw_view(view), brd);
  fill_left_border(rotated90ccw_view(view), brd);
}

// делает свертку входного изображения по строкам с заданным фильтром, помещая его центр 
// в точки, находящиеся к любой из 4-х сторон не ближе, чем brd;
// во входном изображении шагает на 2 пиксела по X
// shift - значение, добавляемое к результирующим писелам, для того чтобы 0 в высокачастотном фильтре выглядел серым
template <typename VI, typename VO>
void convolve_downsample_x(const VI & in, int brdX, int brdY, const VO & out, const std::vector<float> & filter, float shift)
{
  if (brdX < (int)filter.size()/2)
    throw std::runtime_error("too narrow border in convolve_downsample_x");
  if (in.width()/2 - brdX != out.width())
    throw std::runtime_error("half of input image width minus borders at both sides is not equal to output image width");
  if (in.height() - 2 * brdY != out.height())
    throw std::runtime_error("input image height minus borders at both sides is not equal to output image height");

  for (int y = 0; y < out.height(); ++y)
  {
    auto i = in.row_begin(y+brdY) + brdX;
    auto iEnd = in.row_end(y + brdY) - brdX;
    auto o = out.row_begin(y);
    for (; i < iEnd; i += 2, ++o)
    {
      rgb32f_pixel_t sum(shift, shift, shift);
      auto ii = i - (filter.size() - 1) / 2;
      for (auto f : filter)
        sum += f * *ii++;
      *o = sum;
    }
  }
}

// аналог для сертки по столбцам
template <typename VI, typename VO>
void convolve_downsample_y(const VI & in, int brdX, int brdY, const VO & out, const std::vector<float> & filter, float shift)
{
  convolve_downsample_x(transposed_view(in), brdY, brdX, transposed_view(out), filter, shift);
}

// один уровень вейвлет разложения
template <typename VI, typename VO>
void wavelet_transform1(const VI & in, const VO & out, const std::vector<float> & low_pass, const std::vector<float> & hi_pass)
{
  if (in.dimensions() != out.dimensions())
    throw std::runtime_error("input and output images shall have the same dimensions in wavelet_transform");

  // копия входного изображения с дополнительными краями, чтобы фильтры не выходили за пределы
  int brd = std::max(low_pass.size(), hi_pass.size()) / 2;
  rgb32f_image_t extended_in(in.width() + 2 * brd, in.height() + 2 * brd);
  copy_pixels(in, subimage_view(view(extended_in), { brd, brd }, in.dimensions()));
  fill_all_borders(view(extended_in), brd);

  // разложение по X
  rgb32f_image_t filtered_x(in.width() + 2 * brd, in.height() + 2 * brd);
  convolve_downsample_x(view(extended_in), brd, brd,
    subimage_view(view(filtered_x), { brd, brd }, { in.width() / 2, in.height() }), low_pass, 0);
  convolve_downsample_x(view(extended_in), brd, brd,
    subimage_view(view(filtered_x), { brd + in.width() / 2, brd }, { in.width() / 2, in.height() }), hi_pass, 0.5f);
  fill_all_borders(view(filtered_x), brd);

  // разложение по Y
  convolve_downsample_y(view(filtered_x), brd, brd,
    subimage_view(out, { 0, 0 }, { in.width(), in.height() / 2 }), low_pass, 0);
  convolve_downsample_y(view(filtered_x), brd, brd,
    subimage_view(out, { 0, in.height() / 2 }, { in.width(), in.height() / 2 }), hi_pass, 0.5f);
}

// вейвлет разложение с заданным числом уровней
template <typename VI, typename VO>
void wavelet_transform(int levels, const VI & in, const VO & out, const std::vector<float> & low_pass, const std::vector<float> & hi_pass)
{
  if (levels < 1)
    throw std::runtime_error("at least 1 level of transform is required");

  wavelet_transform1(in, out, low_pass, hi_pass);
  if (levels == 1)
    return;

  point2<ptrdiff_t> half_dim(out.width()/2, out.height()/2);
  auto low_freq_out = subimage_view(out, { 0, 0 }, half_dim);

  rgb32f_image_t tmp(half_dim);
  copy_pixels(low_freq_out, view(tmp));
  wavelet_transform(levels-1, view(tmp), low_freq_out, low_pass, hi_pass);
}

// делает свертку входного изображения по строкам с заданным фильтром;
// результат записывается в выходное изобрадение с отсутпом brd от границ;
// в выходном изображении шагает на 2 пиксела по X;
// shift - значение, вычитаемое из входных пикселов, чтобы принимать серый цвет в качестве 0 для высоких частот
template <typename VI, typename VO>
void convolve_upsample_x(const VI & in, int brdX, int brdY, const VO & out, const std::vector<float> & filter, float shift)
{
  if (brdX < (int)filter.size() / 2)
    throw std::runtime_error("too narrow border in convolve_upsample_x");
  if (in.width() != out.width() / 2 - brdX)
    throw std::runtime_error("half of output image width minus borders at both sides is not equal to input image width");
  if (in.height() != out.height() - 2 * brdY)
    throw std::runtime_error("output image height minus borders at both sides is not equal to input image height");

  for (int y = 0; y < in.height(); ++y)
  {
    auto o = out.row_begin(y + brdY) + brdX;
    auto i = in.row_begin(y);
    auto iEnd = in.row_end(y);
    for (; i < iEnd; ++i, o += 2)
    {
      auto inval = *i - rgb32f_pixel_t(shift, shift, shift);
      auto oo = o - (filter.size() - 1) / 2;
      for (auto f : filter)
        *oo++ += f * inval;
    }
  }
}

// аналог для сертки по столбцам
template <typename VI, typename VO>
void convolve_upsample_y(const VI & in, int brdX, int brdY, const VO & out, const std::vector<float> & filter, float shift)
{
  convolve_upsample_x(transposed_view(in), brdY, brdX, transposed_view(out), filter, shift);
}

// один уровень обратного вейвлет разложения
template <typename VI, typename VO>
void inverse_transform1(const VI & in, const VO & out, const std::vector<float> & low_pass, const std::vector<float> & hi_pass)
{
  if (in.dimensions() != out.dimensions())
    throw std::runtime_error("input and output images shall have the same dimensions in inverse_transform");
  auto half_width = out.width() / 2;
  auto half_height = out.height() / 2;

  // обратное преобразование по Y
  int brd = std::max(low_pass.size(), hi_pass.size()) / 2;
  rgb32f_image_t inverted_y(in.width(), in.height() + 2 * brd);
  fill_pixels(view(inverted_y), rgb32f_pixel_t(0, 0, 0));
  convolve_upsample_y(subimage_view(in, { 0, 0 }, { in.width(), half_height }),
    0, brd, view(inverted_y), low_pass, 0);
  convolve_upsample_y(subimage_view(in, { 0, half_height }, { in.width(), half_height }),
    0, brd, view(inverted_y), hi_pass, 0.5f);

  // обратное преобразование по X
  rgb32f_image_t inverted_x(in.width() + 2 * brd, in.height());
  fill_pixels(view(inverted_x), rgb32f_pixel_t(0, 0, 0));
  convolve_upsample_x(subimage_view(const_view(inverted_y), { 0, brd }, { half_width, in.height() }),
    brd, 0, view(inverted_x), low_pass, 0);
  convolve_upsample_x(subimage_view(const_view(inverted_y), { half_width, brd }, { half_width, in.height() }),
    brd, 0, view(inverted_x), hi_pass, 0.5f);

  // копирование результата без границ
  copy_pixels(subimage_view(const_view(inverted_x), { brd, 0 }, out.dimensions()), out);
}

// обратное вейвлет разложение с заданным числом уровней
template <typename VI, typename VO>
void inverse_transform(int levels, const VI & in, const VO & out, const std::vector<float> & low_pass, const std::vector<float> & hi_pass)
{
  if (levels < 1)
    throw std::runtime_error("at least 1 level of transform is required");

  if (levels == 1)
  {
    inverse_transform1(in, out, low_pass, hi_pass);
    return;
  }

  point2<ptrdiff_t> half_dim(out.width() / 2, out.height() / 2);

  rgb32f_image_t tmp(in.dimensions());
  copy_pixels(in, view(tmp));
  inverse_transform(levels - 1, subimage_view(in, { 0, 0 }, half_dim), subimage_view(view(tmp), { 0, 0 }, half_dim), low_pass, hi_pass);

  inverse_transform1(view(tmp), out, low_pass, hi_pass);
}

// демонстрация прямого и обратного вейвлет преобразований с записью результатов в файлы
template <typename V>
void demo_transform(const V & img, const std::string & name, 
  const std::vector<float> & low_pass_synthesis, const std::vector<float> & hi_pass_synthesis)
{
  // фильтры отличаются на множетель 2, чтобы низкие частоты прямого преобразования выглядели как усреднение
  auto low_pass_analysis = low_pass_synthesis;
  for (auto & v : low_pass_analysis)
    v /= 2;
  auto hi_pass_analysis = hi_pass_synthesis;
  for (auto & v : hi_pass_analysis)
    v /= 2;

  rgb32f_image_t transformed(img.dimensions());
  wavelet_transform(3, img, view(transformed), low_pass_analysis, hi_pass_analysis);
  png_write_float_view(("transformed-" + name + ".png").c_str(), const_view(transformed));

  rgb32f_image_t restored(img.dimensions());
  inverse_transform(3, const_view(transformed), view(restored), low_pass_synthesis, hi_pass_synthesis);
  png_write_float_view(("restored-" + name + ".png").c_str(), const_view(restored));
}

void main()
{
  // считываем объект
  rgb32f_image_t img;
  png_read_float_image("lena.png", img);

  // https://en.wikipedia.org/wiki/Daubechies_wavelet
  demo_transform(const_view(img), "D2", { 1, 1 }, { 1, -1 }); //Haar

  demo_transform(const_view(img), "D4", 
    { 0.6830127f, 1.1830127f, 0.3169873f, -0.1830127f }, 
    { -0.1830127f, -0.3169873f, 1.1830127f, -0.6830127f });

  demo_transform(const_view(img), "D6",
    { 0.47046721f, 1.14111692f, 0.650365f, -0.19093442f, -0.12083221f, 0.0498175f },
    { 0.0498175f, 0.12083221f, -0.19093442f, -0.650365f, 1.14111692f, -0.47046721f });

  // LeGall - something wrong with hi-freqs restoration
/*
  wavelet_transform1(const_view(img), view(transformed), { -0.125f, 0.25f, 0.75f, 0.25f, -0.125f }, { -0.5f, 1.0f, -0.5f });
  png_write_float_view("legall.png", const_view(transformed));

  inverse_transform1(const_view(transformed), view(restored), { 0.5f, 1.0f, 0.5f }, { -0.125f, -0.25f, 0.75f, -0.25f, -0.125f });
  png_write_float_view("restored-legall.png", const_view(restored));*/
}
