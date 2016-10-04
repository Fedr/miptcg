#pragma once

// арифметические операции с пикселами (отдельно в каждом канале)

inline rgb32f_pixel_t & operator += (rgb32f_pixel_t & a, const rgb32f_pixel_t & b)
{
  get_color(a, red_t()) += get_color(b, red_t());
  get_color(a, green_t()) += get_color(b, green_t());
  get_color(a, blue_t()) += get_color(b, blue_t());
  return a;
}

inline rgb32f_pixel_t operator + (rgb32f_pixel_t a, const rgb32f_pixel_t & b)
{
  return rgb32f_pixel_t(a) += b;
}

inline rgb32f_pixel_t & operator -= (rgb32f_pixel_t & a, const rgb32f_pixel_t & b)
{
  get_color(a, red_t()) -= get_color(b, red_t());
  get_color(a, green_t()) -= get_color(b, green_t());
  get_color(a, blue_t()) -= get_color(b, blue_t());
  return a;
}

inline rgb32f_pixel_t operator - (rgb32f_pixel_t a, const rgb32f_pixel_t & b)
{
  return rgb32f_pixel_t(a) -= b;
}

inline rgb32f_pixel_t operator * (float a, rgb32f_pixel_t b)
{
  get_color(b, red_t()) *= a;
  get_color(b, green_t()) *= a;
  get_color(b, blue_t()) *= a;
  return b;
}

inline float operator * (const rgb32f_pixel_t & a, const rgb32f_pixel_t & b)
{
  return
    get_color(a, red_t()) * get_color(b, red_t()) +
    get_color(a, green_t()) * get_color(b, green_t()) +
    get_color(a, blue_t()) * get_color(b, blue_t());
}
