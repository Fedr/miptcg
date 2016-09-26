#pragma once

// запись картинки с пикселами в представлении с плавающей точкой в PNG-формат
// значения меньше 0 записываются как минимальная интесивность 0
// значения больше 1 записываются как максимальная интесивность 255
template <typename View>
inline void png_write_float_view(const char* filename, const View& view)
{
  png_write_view(filename, color_converted_view<rgb8_pixel_t>(view,
    [](const auto & src, auto & dst) { static_for_each(src, dst, [](auto f, auto & i) { i = (f < 0) ? 0 : (f > 1 ? 255 : int(f*255.5f)); }); }));
}

// загружает изображение как rgb32f_image_t
inline void png_read_float_image(const char* filename, rgb32f_image_t& img)
{
  rgb8_image_t img8;
  png_read_image(filename, img8);
  img.recreate(img8.dimensions());
  copy_pixels(color_converted_view<rgb32f_pixel_t>(const_view(img8)), view(img));
}
