//по материалам статьи
//https://www.cs.jhu.edu/~misha/Fall07/Papers/Perez03.pdf
//чтобы собрать программу, надо слинковать с libpng.lib и zlib.lib

#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_io.hpp>
using namespace boost::gil;

// копирование пикселов маски из from в to
template <typename M, typename V, typename VT>
void clone(const M & mask, const V & from, const VT & to)
{
  if (from.dimensions() != mask.dimensions() || from.dimensions() != to.dimensions())
    throw std::runtime_error("image dimensions shall be equal");
  auto ito = to.begin();
  auto imask = mask.begin();
  for (const auto & p : from)
  {
    if (*imask)
      *ito = p;
    ++ito;
    ++imask;
  }
}

// арифметические операции с пикселами (отдельно в каждом канале)
inline rgb32f_pixel_t operator + (rgb32f_pixel_t a, const rgb32f_pixel_t & b)
{
  get_color(a, red_t()) += get_color(b, red_t());
  get_color(a, green_t()) += get_color(b, green_t());
  get_color(a, blue_t()) += get_color(b, blue_t());
  return a;
}

inline rgb32f_pixel_t operator - (rgb32f_pixel_t a, const rgb32f_pixel_t & b)
{
  get_color(a, red_t()) -= get_color(b, red_t());
  get_color(a, green_t()) -= get_color(b, green_t());
  get_color(a, blue_t()) -= get_color(b, blue_t());
  return a;
}

inline float absmax(float a, float b)
{
  return fabs(a) >= fabs(b) ? a : b;
}

inline rgb32f_pixel_t absmax(rgb32f_pixel_t a, const rgb32f_pixel_t & b)
{
  get_color(a, red_t()) = absmax(get_color(a, red_t()), get_color(b, red_t()));
  get_color(a, green_t()) = absmax(get_color(a, green_t()), get_color(b, green_t()));
  get_color(a, blue_t()) = absmax(get_color(a, blue_t()), get_color(b, blue_t()));
  return a;
}

inline rgb32f_pixel_t operator * (float a, rgb32f_pixel_t b)
{
  get_color(b, red_t()) *= a;
  get_color(b, green_t()) *= a;
  get_color(b, blue_t()) *= a;
  return b;
}

// запоминает смещения до четырех соседних пикселов
template <typename L>
struct cross_locations
{
  cross_locations(const L & loc)
    : w(loc.cache_location(-1, 0))
    , n(loc.cache_location(0, 1))
    , s(loc.cache_location(0, -1))
    , e(loc.cache_location(1, 0))
  {
  }
  typename L::cached_location_t w, n, s, e;
};

template <typename L>
inline cross_locations<L> cross(const L & loc)
{
  return { loc };
}

// одна итерация решения задачи Пуассона методом Гаусса-Зейделя:
// https://ru.wikipedia.org/wiki/Метод_Гаусса_—_Зейделя_решения_системы_линейных_уравнений
// mask - только точки, выбранные маской, будут меняться
// sol - исходное и результирующее приближения решения
// rhs - правая часть уравнения
template <typename M, typename S, typename R>
void poisson1(const M & mask, const S & sol, const R & rhs)
{
  if (sol.dimensions() != mask.dimensions() || sol.dimensions() != rhs.dimensions())
    throw std::runtime_error("image dimensions shall be equal");

  auto sol_cross = cross(sol.xy_at(0, 0));

  for (int y = 1; y + 1 < sol.height(); ++y)
  {
    auto sol_loc = sol.xy_at(1, y);
    auto imask = std::next(mask.row_begin(y));
    auto irhs = std::next(rhs.row_begin(y));
    for (int x = 1; x + 1 < sol.width(); ++x, ++imask, ++irhs, ++sol_loc.x())
    {
      if (!*imask)
        continue;
      *sol_loc = 0.25f * (
        sol_loc[sol_cross.w] + sol_loc[sol_cross.n] + sol_loc[sol_cross.s] + sol_loc[sol_cross.e] - *irhs);
    }
  }
}

// n итераций методом Гаусса-Зейделя
template <typename M, typename S, typename R>
void poisson(int n, const M & mask, const S & sol, const R & rhs)
{
  for (int i = 0; i < n; ++i)
    poisson1(mask, sol, rhs);
}

// вычисляет лапласиан данного изображения в каждой точке маски
template <typename M, typename V, typename L>
void get_laplacian(const M & mask, const V & img, const L & laplacian)
{
  if (img.dimensions() != mask.dimensions() || img.dimensions() != laplacian.dimensions())
    throw std::runtime_error("image dimensions shall be equal");

  auto img_cross = cross(img.xy_at(0, 0));

  for (int y = 1; y + 1 < img.height(); ++y)
  {
    auto img_loc = img.xy_at(1, y);
    auto imask = std::next(mask.row_begin(y));
    auto il = std::next(laplacian.row_begin(y));
    for (int x = 1; x + 1 < img.width(); ++x, ++imask, ++il, ++img_loc.x())
    {
      if (!*imask)
        continue;
      *il = img_loc[img_cross.w] + img_loc[img_cross.n] + img_loc[img_cross.s] + img_loc[img_cross.e] - 4 * *img_loc;
    }
  }
}

// вычисляет "лапласиан", используя максимальные по модулю разности из двух изображений
template <typename M, typename V1, typename V2, typename L>
void get_absmax_laplacian(const M & mask, const V1 & img1, const V2 & img2, const L & laplacian)
{
  if (img1.dimensions() != mask.dimensions() || img1.dimensions() != laplacian.dimensions() || img2.dimensions() != laplacian.dimensions())
    throw std::runtime_error("image dimensions shall be equal");

  auto img1_cross = cross(img1.xy_at(0, 0));
  auto img2_cross = cross(img2.xy_at(0, 0));

  for (int y = 1; y + 1 < img1.height(); ++y)
  {
    auto img1_loc = img1.xy_at(1, y);
    auto img2_loc = img2.xy_at(1, y);
    auto imask = std::next(mask.row_begin(y));
    auto il = std::next(laplacian.row_begin(y));
    for (int x = 1; x + 1 < img1.width(); ++x, ++imask, ++il, ++img1_loc.x(), ++img2_loc.x())
    {
      if (!*imask)
        continue;
      *il = 
        absmax(img1_loc[img1_cross.w] - *img1_loc, img2_loc[img2_cross.w] - *img2_loc) +
        absmax(img1_loc[img1_cross.n] - *img1_loc, img2_loc[img2_cross.n] - *img2_loc) +
        absmax(img1_loc[img1_cross.s] - *img1_loc, img2_loc[img2_cross.s] - *img2_loc) +
        absmax(img1_loc[img1_cross.e] - *img1_loc, img2_loc[img2_cross.e] - *img2_loc);
    }
  }
}

// помечаем пикселы маски рядом с правой границей числом 2
template <typename M>
void mark_xpos(const M & mask)
{
  for (int y = 0; y < mask.height(); ++y)
  { 
    auto i = mask.row_begin(y), iEnd = mask.row_end(y);
    bool prev = false;
    for (; i != iEnd; ++i)
    {
      bool curr = *i != 0;
      if (!prev && curr)
        *i = 2;
      prev = curr;
    }
  }
}

// выедание пикселов маски со всех четырех сторон
template <typename M>
void erode(const M & mask)
{
  mark_xpos(mask);
  mark_xpos(rotated180_view(mask));
  mark_xpos(rotated90cw_view(mask));
  mark_xpos(rotated90ccw_view(mask));
  //замена 2 на 0
  for (auto & p : mask) 
    if (p == 2) 
      p = 0;
}

// запись картинки с пикселами в представлении с плавающей точкой в PNG-формат
// значения меньше 0 записываются как минимальная интесивность 0
// значения больше 1 записываются как максимальная интесивность 255
template <typename View>
inline void png_write_float_view(const char* filename, const View& view)
{
  png_write_view(filename, color_converted_view<rgb8_pixel_t>(view,
    [](const auto & src, auto & dst) { static_for_each(src, dst, [](auto f, auto & i) { i = (f < 0) ? 0 : (f > 1 ? 255 : int(f*255.5f)); }); }));
}

// функтор, который возвращает (0,0,0) для любой точки
struct zero
{
  using argument_type = point2<ptrdiff_t>;
  using value_type = rgb32f_pixel_t;
  using reference  = rgb32f_pixel_t;
  using const_t = zero;
  rgb32f_pixel_t operator()(const point2<ptrdiff_t>& p) const { return{ 0, 0, 0 }; }
};

void main()
{
  // считываем объект и переводи его в представление с плавающей точкой
  rgb8_image_t fore;
  png_read_image("fore.png", fore);
  rgb32f_image_t foref(fore.dimensions());
  copy_pixels(color_converted_view<rgb32f_pixel_t>(const_view(fore)), view(foref));

  // считываем фон и переводи его в представление с плавающей точкой
  rgb8_image_t back;
  png_read_image("back.png", back);
  rgb32f_image_t backf(back.dimensions());
  copy_pixels(color_converted_view<rgb32f_pixel_t>(const_view(back)), view(backf));

  //определям маску как не-черные точки объекта
  gray8_image_t mask(fore.dimensions());
  copy_pixels(color_converted_view<gray8_pixel_t>(const_view(fore), 
    [](const auto & src, auto & dst) { dst = (get_color(src, red_t()) != 0 || get_color(src, green_t()) != 0 || get_color(src, blue_t())) != 0 ? 1 : 0; }), view(mask));
  //уменьшаем маску на 1 пиксель со всех сторон, чтобы можно было посчитать градиент объекта во всех точках маски
  erode(view(mask));

  //простое клонирование объекта
  rgb32f_image_t clonef = backf;
  clone(const_view(mask), const_view(foref), view(clonef));
  png_write_float_view("clone.png", const_view(clonef));

  // заполнение дырки в фоне методом Лапласа (объект игнорируется)
  rgb32f_image_t laplacef = backf;
  using zero_locator = virtual_2d_locator<zero, false>;
  image_view<zero_locator> zero_rhs(backf.dimensions(), zero_locator());
  poisson(300, const_view(mask), view(laplacef), zero_rhs);
  png_write_float_view("laplace.png", const_view(laplacef));

  // заполнение дырки в фоне, копируя градиент из объекта
  rgb32f_image_t importf = backf;
  rgb32f_image_t fore_laplacian(foref.dimensions());
  get_laplacian(const_view(mask), const_view(foref), view(fore_laplacian));
  poisson(300, const_view(mask), view(importf), view(fore_laplacian));
  png_write_float_view("import.png", const_view(importf));

  // заполнение дырки в фоне, используя максимальный градиент из объекта или фона
  rgb32f_image_t mixedf = backf;
  rgb32f_image_t max_fore_back_laplacian(foref.dimensions());
  get_absmax_laplacian(const_view(mask), const_view(foref), const_view(backf), view(max_fore_back_laplacian));
  poisson(300, const_view(mask), view(mixedf), view(max_fore_back_laplacian));
  png_write_float_view("mixed.png", const_view(mixedf));
}
