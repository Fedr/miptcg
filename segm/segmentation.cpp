// (c) Fedor Chelnokov, fchel@mail.ru

#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/jpeg_io.hpp>

#define DIR "C:\\graphics\\images\\"
#define FILENAME DIR "42049"
const float MU = 5;
const float NU = 100;
const float GAMMA2 = 0.3f * 0.3f;
const float TETHA_E = 10; //more white
const float TETHA_D = 10; //more black

// define images and view, those pixels are floats (unlike gray32f_pixel_t, which float limited in the range [0,1])
namespace boost { namespace gil {
  typedef float bits32fu;
  GIL_DEFINE_BASE_TYPEDEFS(32fu,gray)
} }

using namespace boost::gil;

// this type is implicitly convertible to/from gray32f_pixel_t and to/from float
typedef channel_type<gray32f_view_t>::type pixel_float_t;

// creates view that applies given functor to the base view
template <typename ViewType, typename Functor>
inline typename ViewType::template add_deref<Functor>::type function_view(const ViewType & view, const Functor & f)
{
  return ViewType::template add_deref<Functor>::make(view, f);
}

// compute prior probability be the formula below Fig.6
void find_prior_probability(const gray8c_view_t & pic, const gray32f_view_t & prob)
{
  transform_pixels(pic, prob,
    [](unsigned char c) -> float 
    { 
      if (c == 0)
        return 0;
      if (c == 255)
        return 1;
      float p1 = (1.0f / 255) * c;
      float p0 = 1 - p1;
      float t = log(p1 / p0);
      return 1 / (1 + exp(-t / MU)); 
    }
  );
};

// computes square of the gradient
inline float sqr_diff(float a, float b)
{
  float d = a - b;
  return d * d;
}

// if given value is less than stored in cell, then updates cell and increases the counter
inline void updateDistance(float value, gray32fu_pixel_t & cell, int & changes)
{
  if (value < cell)
  {
    cell = value;
    ++changes;
  }
}

// makes one pass from left to right from top to bottom, updating the distance transform
template <typename PicView, typename DView>
int improve_ggdt_forward(const PicView & pic, const DView & d)
{
  int changes = 0;

  auto pic_loc = pic.xy_at(1, 0);
  const auto pic_w = pic_loc.cache_location(-1,0);
  const auto pic_nw = pic_loc.cache_location(-1,-1);
  const auto pic_n = pic_loc.cache_location(0,-1);
  const auto pic_ne = pic_loc.cache_location(1,-1);

  auto d_loc = d.xy_at(1, 0);
  const auto d_w = d_loc.cache_location(-1,0);
  const auto d_nw = d_loc.cache_location(-1,-1);
  const auto d_n = d_loc.cache_location(0,-1);
  const auto d_ne = d_loc.cache_location(1,-1);

  //process first row
  for (int x = 1; x < pic.width(); ++x, ++pic_loc.x(), ++d_loc.x())
  {
    updateDistance(d_loc[d_w] + sqrt(1 + GAMMA2 * sqr_diff(pic_loc[pic_w], *pic_loc)),
      *d_loc, changes);
  }
  
  //process other rows
  for (int y = 1; y < pic.height(); ++y)
  {
    pic_loc = pic.xy_at(0, y);
    d_loc = d.xy_at(0, y);

    //process first column
    updateDistance(
      std::min(
        d_loc[d_n] + sqrt(1 + GAMMA2 * sqr_diff(pic_loc[pic_n], *pic_loc)),
        d_loc[d_ne] + sqrt(2 + GAMMA2 * sqr_diff(pic_loc[pic_ne], *pic_loc))),
      *d_loc, changes);

    //process other columns
    ++pic_loc.x(), ++d_loc.x();
    for (int x = 1; x < pic.width(); ++x, ++pic_loc.x(), ++d_loc.x())
    {
      updateDistance(
        std::min(std::min(std::min(
          d_loc[d_w] + sqrt(1 + GAMMA2 * sqr_diff(pic_loc[pic_w], *pic_loc)),
          d_loc[d_nw] + sqrt(2 + GAMMA2 * sqr_diff(pic_loc[pic_nw], *pic_loc))),
          d_loc[d_n] + sqrt(1 + GAMMA2 * sqr_diff(pic_loc[pic_n], *pic_loc))),
          d_loc[d_ne] + sqrt(2 + GAMMA2 * sqr_diff(pic_loc[pic_ne], *pic_loc))),
        *d_loc, changes);
    }
  }

  return changes;
}

// makes the pass from left to right from top to bottom, updating the distance transform,
// and the pass in the opposite direction
int improve_ggdt(const gray8c_view_t & pic, const gray32fu_view_t & d)
{
  return
    improve_ggdt_forward(pic, d) +
    improve_ggdt_forward(rotated180_view(pic), rotated180_view(d));
}

// computes generalized geodesic distance 
template <typename ProbView>
void find_ggdt(const gray8c_view_t & pic, const ProbView & prob, const gray32fu_view_t & d)
{
  //scale seed mask
  transform_pixels(prob, d, [](pixel_float_t v) -> float { return NU * v; } );

  const int MAX_ITERS = 10;
  for (int i = 0; i < MAX_ITERS; ++i)
  {
    int changes = improve_ggdt(pic, d);
    if (changes == 0)
      return;
  }
}

// gray32f_pixel_t -> gray32f_pixel_t: y = 1 - x
struct completer : deref_base<completer, gray32f_pixel_t, gray32f_pixel_t, const gray32f_pixel_t&, gray32f_pixel_t, gray32f_pixel_t, false> 
{
  gray32f_pixel_t operator()(pixel_float_t x) const 
    { return gray32f_pixel_t(1 - x); }
};

// computes signed generalized geodesic distance
void find_ds(const gray8c_view_t & pic, const gray32fc_view_t & prob, const gray32fu_view_t & ds)
{
  find_ggdt(pic, prob, ds);

  gray32fu_image_t d(ds.dimensions());
  find_ggdt(pic, function_view(prob, completer()), view(d));

  //ds -= d;
  transform_pixels(ds, const_view(d), ds,
    [](float a, float b) -> float { return a - b; }
  );
}

// gray32fu_pixel_t -> gray32f_pixel_t: computes y = (x > t) ? a : b
class discretizor : public deref_base<discretizor, gray32f_pixel_t, gray32f_pixel_t, const gray32f_pixel_t&, float, gray32f_pixel_t, false> 
{
  float t, a, b;
public:
  discretizor(float t, float a, float b) : t(t), a(a), b(b) { }
  gray32f_pixel_t operator()(float x) const 
    { return gray32f_pixel_t(x > t ? a : b); }
};

// computes symmetric signed distance
template <typename MView>
void find_dss(const gray8c_view_t & pic, const MView & Me, const MView & notMd, const gray32fu_view_t & dss)
{
  find_ggdt(pic, Me, dss);

  gray32fu_image_t d(dss.dimensions());
  find_ggdt(pic, notMd, view(d));

  //dss -= d + TETHA_D - TETHA_E;
  transform_pixels(dss, const_view(d), dss,
    [](float a, float b) -> float { return a - b + TETHA_D - TETHA_E; }
  );
}

// writes view containing positive and negative values such as 0 becomes 127-gray value, and choosing appropriate scale
void jpeg_normalized_write_view(const char * filename, const gray32fuc_view_t & view)
{
  float maxd = *std::max_element(view.begin(), view.end());
  float mind = *std::min_element(view.begin(), view.end());
  float fact = std::max(maxd, -mind);
  jpeg_write_view(filename, color_converted_view<gray8_pixel_t>(view,
    [fact](float in, gray8_pixel_t & out) { out = unsigned char(((in + fact) * 255 / (2 * fact)) + 0.5f); } ));
}

void main()
{
  gray8_image_t image;
  jpeg_read_image(FILENAME ".jpg", image);
  auto dim = image.dimensions();

  gray32f_image_t prob(dim);
  find_prior_probability(const_view(image), view(prob));
  jpeg_write_view(FILENAME "-1-prob.jpg", color_converted_view<gray8_pixel_t>(view(prob)));

  gray32fu_image_t ds(dim);
  find_ds(const_view(image), const_view(prob), view(ds));

  jpeg_normalized_write_view(FILENAME "-2-ds.jpg", const_view(ds));

  auto Me = function_view(const_view(ds), discretizor(-TETHA_E, 1, 0));
  auto notMd = function_view(const_view(ds), discretizor(TETHA_D, 0, 1));

  jpeg_write_view(FILENAME "-3-Md-not.jpg", color_converted_view<gray8_pixel_t>(notMd));
  jpeg_write_view(FILENAME "-4-Me.jpg", color_converted_view<gray8_pixel_t>(Me));

  gray32fu_image_t dss(dim);
  find_dss(const_view(image), Me, notMd, view(dss));
  jpeg_normalized_write_view(FILENAME "-5-dss.jpg", const_view(dss));

  auto segm = function_view(const_view(dss), discretizor(0, 1, 0));
  jpeg_write_view(FILENAME "segm.jpg", color_converted_view<gray8_pixel_t>(segm));
}
