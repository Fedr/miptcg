#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/jpeg_io.hpp>
#include <boost/gil/extension/io/png_io.hpp>
using namespace boost::gil;

#include <vector>
#include <algorithm>

using pix_values = std::vector<float>;
using pix_indices = std::vector<int>;

pix_indices get_asc_order_indices(const pix_values & i_img)
{
  pix_indices res(i_img.size());
  for (size_t i = 0; i < i_img.size(); ++i)
    res[i] = i;
  std::sort(res.begin(), res.end(), [&i_img](int a, int b) { return i_img[a] < i_img[b]; });
  return res;
}

void copy_1d_hist(pix_values & io_img, const pix_values & i_target)
{
  if (io_img.size() != i_target.size())
    throw std::runtime_error("images of differen sizes are not supported");
  auto ind0 = get_asc_order_indices(io_img);
  auto indT = get_asc_order_indices(i_target);

  const float INERTIA = 0.75f;
  for (size_t i = 0; i < io_img.size(); ++i)
  { 
    auto & v = io_img[ind0[i]];
    v = INERTIA * v + (1 - INERTIA) * i_target[indT[i]];
  }
}

struct color_dir
{
  float r;
  float g;
  float b;

  color_dir(float ir, float ig, float ib)
  {
    auto rlen = 1 / sqrt(ir*ir + ig*ig + ib*ib);
    r = ir * rlen;
    g = ig * rlen;
    b = ib * rlen;
  }

  template <typename P>
  float get_proj(const P & p) const
  {
    return r * get_color(p, red_t()) + g * get_color(p, green_t()) + b * get_color(p, blue_t());
  }

  template <typename P>
  P set_proj(const P & p, float val) const
  {
    P res = p;
    val -= get_proj(p);
    get_color(res, red_t()) += r * val;
    get_color(res, green_t()) += g * val;
    get_color(res, blue_t()) += b * val;
    return res;
  }
};

template <typename V>
pix_values get_pix_values_in_color_direction(const V & view, const color_dir & dir)
{
  pix_values res;
  res.reserve(view.width() * view.height());

  for (const auto & p : view)
  {
    res.push_back(dir.get_proj(p));
  }

  assert(res.size() == view.width() * view.height());
  return res;
}

template <typename V>
void set_pix_values_in_color_direction(const V & view, const pix_values & vals, const color_dir & dir)
{
  if (vals.size() != view.width() * view.height())
    throw std::runtime_error("wrong number of values");

  auto it = vals.begin();
  for (auto & p : view)
  {
    p = dir.set_proj(p, *it++);
  }
}

template <typename V, typename VT>
void copy_hist_in_dir(const V & img, const VT & target, const color_dir & dir)
{
  auto vals = get_pix_values_in_color_direction(img, dir);
  const auto valsTarget = get_pix_values_in_color_direction(target, dir);
  copy_1d_hist(vals, valsTarget);
  set_pix_values_in_color_direction(img, vals, dir);
}

void main()
{
  rgb8_image_t a;
  jpeg_read_image("a.jpg", a);

  rgb8_image_t b;
  jpeg_read_image("b.jpg", b);

  rgb32f_image_t af(a.dimensions());
  copy_pixels(color_converted_view<rgb32f_pixel_t>(const_view(a)), view(af));

  rgb32f_image_t bf(b.dimensions());
  copy_pixels(color_converted_view<rgb32f_pixel_t>(const_view(b)), view(bf));

  copy_hist_in_dir(view(af), view(bf), { 1,0,0 });
  copy_hist_in_dir(view(af), view(bf), { 0,1,0 });
  copy_hist_in_dir(view(af), view(bf), { 0,0,1 });

/*  std::srand(0);
  for (int i = 0; i < 300; ++i)
  {
    color_dir dir(float(std::rand()) - RAND_MAX / 2.0f, float(std::rand()) - RAND_MAX / 2.0f, float(std::rand()) - RAND_MAX / 2.0f);
    copy_hist_in_dir(view(af), view(bf), dir);
  }*/

  png_write_view("a-res.png", color_converted_view<rgb8_pixel_t>(const_view(af)));
  jpeg_write_view("a-res.jpg", color_converted_view<rgb8_pixel_t>(const_view(af)));
}
