#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;
    this->sample_unit = 1.0 / sqrt(this->sample_rate);
    this->sample_sqrt = int(sqrt(this->sample_rate));

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(float x, float y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    if (x < 0 || x >= width) return;
    if (y < 0 || y >= height) return;
    
    int j = floor(x), i = floor(y);
    int nx = (x-j) * sample_sqrt;
    int ny = (y-i) * sample_sqrt;
    sample_buffer[(i*sample_sqrt+ny) * width * sample_sqrt + j*sample_sqrt + nx] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    // check bounds
    int sx = floor(x), sy = floor(y);
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;
    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  bool RasterizerImp::p_in_triangle(float x, float y, float x0, float y0, float x1, float y1, float x2, float y2) {
    CGL::Vector3D e0 = CGL::Vector3D(x1-x0, y1-y0, 0); // vector from point 0
    CGL::Vector3D e1 = CGL::Vector3D(x2-x1, y2-y1, 0); // vector from point 1
    CGL::Vector3D e2 = CGL::Vector3D(x0-x2, y0-y2, 0); // vector from point 2

    CGL::Vector3D v0 = CGL::cross(e0, CGL::Vector3D(x-x0, y-y0, 0));
    CGL::Vector3D v1 = CGL::cross(e1, CGL::Vector3D(x-x1, y-y1, 0));
    CGL::Vector3D v2 = CGL::cross(e2, CGL::Vector3D(x-x2, y-y2, 0));

    return (v0[2]*v1[2] >= 0 && v1[2]*v2[2] >= 0 && v2[2]*v0[2] >= 0);
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // TODO: Task 2: Update to implement super-sampled rasterization
    float lBound, rBound, tBound, bBound;
    lBound = min(min(x0, x1), x2), rBound = max(max(x0, x1), x2);
    bBound = min(min(y0, y1), y2), tBound = max(max(y0, y1), y2);
    for (int i = floor(bBound); i <= floor(tBound); i++)
      for (int j = floor(lBound); j <= floor(rBound); ++j)
        for (float y = i; y < i+1; y += sample_unit)
          for (float x = j; x < j+1; x += sample_unit)
            if (p_in_triangle(x+0.5*sample_unit, y+0.5*sample_unit, x0, y0, x1, y1, x2, y2))
              fill_pixel(x, y, color);
    
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    float lBound, rBound, tBound, bBound;
    lBound = min(min(x0, x1), x2), rBound = max(max(x0, x1), x2);
    bBound = min(min(y0, y1), y2), tBound = max(max(y0, y1), y2);
    for (int i = floor(bBound); i <= floor(tBound); i++)
      for (int j = floor(lBound); j <= floor(rBound); ++j) {
        Color color;
        for (float y = i; y < i+1; y += sample_unit)
          for (float x = j; x < j+1; x += sample_unit)
            if (p_in_triangle(x+0.5*sample_unit, y+0.5*sample_unit, x0, y0, x1, y1, x2, y2)) {
              float alpha = ( -(x-x1)*(y2-y1) + (y-y1)*(x2-x1) )/
                            ( -(x0-x1)*(y2-y1) + (y0-y1)*(x2-x1) );
              float beta = ( -(x-x2)*(y0-y2) + (y-y2)*(x0-x2) )/
                           ( -(x1-x2)*(y0-y2) + (y1-y2)*(x0-x2) );
              float gamma = 1 - alpha - beta;
              color = alpha*c0 + beta*c1 + gamma*c2;
              fill_pixel(x, y, color);
            }
      }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
    float lBound, rBound, tBound, bBound;
    lBound = min(min(x0, x1), x2), rBound = max(max(x0, x1), x2);
    bBound = min(min(y0, y1), y2), tBound = max(max(y0, y1), y2);
    for (int i = floor(bBound); i <= floor(tBound); i++)
      for (int j = floor(lBound); j <= floor(rBound); ++j) {
        for (float y = i; y < i+1; y += sample_unit)
          for (float x = j; x < j+1; x += sample_unit) {
            if (!p_in_triangle(x+0.5*sample_unit, y+0.5*sample_unit, x0, y0, x1, y1, x2, y2))
              continue;

            float alpha = ( -(x+0.5*sample_unit-x1)*(y2-y1) + (y+0.5*sample_unit-y1)*(x2-x1) )/
                          ( -(x0-x1)*(y2-y1) + (y0-y1)*(x2-x1) );
            float beta = ( -(x+0.5*sample_unit-x2)*(y0-y2) + (y+0.5*sample_unit-y2)*(x0-x2) )/
                          ( -(x1-x2)*(y0-y2) + (y1-y2)*(x0-x2) );
            float gamma = 1 - alpha - beta;

            float alpha_dx = ( -(x+1.5*sample_unit-x1)*(y2-y1) + (y+0.5*sample_unit-y1)*(x2-x1) )/
                          ( -(x0-x1)*(y2-y1) + (y0-y1)*(x2-x1) );
            float beta_dx = ( -(x+1.5*sample_unit-x2)*(y0-y2) + (y+0.5*sample_unit-y2)*(x0-x2) )/
                          ( -(x1-x2)*(y0-y2) + (y1-y2)*(x0-x2) );
            float gamma_dx = 1 - alpha_dx - beta_dx;

            float alpha_dy = ( -(x+0.5*sample_unit-x1)*(y2-y1) + (y+1.5*sample_unit-y1)*(x2-x1) )/
                          ( -(x0-x1)*(y2-y1) + (y0-y1)*(x2-x1) );
            float beta_dy = ( -(x+0.5*sample_unit-x2)*(y0-y2) + (y+1.5*sample_unit-y2)*(x0-x2) )/
                          ( -(x1-x2)*(y0-y2) + (y1-y2)*(x0-x2) );
            float gamma_dy = 1 - alpha_dy - beta_dy;
            
            SampleParams sp;
            sp.lsm = lsm;
            sp.psm = psm;
            sp.p_uv = alpha*Vector2D(u0, v0) + beta*Vector2D(u1, v1) + gamma*Vector2D(u2, v2);
            sp.p_dx_uv = alpha_dx*Vector2D(u0, v0) + beta_dx*Vector2D(u1, v1) + gamma_dx*Vector2D(u2, v2);
            sp.p_dy_uv = alpha_dy*Vector2D(u0, v0) + beta_dy*Vector2D(u1, v1) + gamma_dy*Vector2D(u2, v2);
            
            fill_pixel(x, y, tex.sample(sp));
          }
        
      }
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support
    this->sample_rate = rate;
    this->sample_unit = 1.0 / sqrt(sample_rate);
    this->sample_sqrt = int(sqrt(sample_rate));
    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }

  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col;
        for (int i = 0; i < sample_sqrt; i++)
          for (int j = 0; j < sample_sqrt; j++) {
            col += sample_buffer[(y*sample_sqrt + i) * width * sample_sqrt + x*sample_sqrt + j];
          }
        col *= 1.0/sample_rate;
        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }
  }

  Rasterizer::~Rasterizer() { }


}// CGL

