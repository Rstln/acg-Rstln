#include <filesystem>
// #include <experimental/filesystem> // uncomment here if the <filesystem> cannot be included above
//
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "Eigen/Core"
//
#include "parse_svg.h"

#ifndef eps
#define eps 1e-2
#endif

/***
 * signed area of a triangle connecting points (p0, p1, p2) in counter-clockwise order.
 * @param p0 1st point xy-coordinate
 * @param p1 2nd point xy-coordinate
 * @param p2 3rd point xy-coordinate
 * @return signed area (float)
 */
float area(
    const Eigen::Vector2f &p0,
    const Eigen::Vector2f &p1,
    const Eigen::Vector2f &p2) {
  const auto v01 = p1 - p0;
  const auto v02 = p2 - p0;
  // return 0.5f * (v01[0] * v02[1] - v01[1] * v02[0]); // right handed coordinate
  return 0.5f * (v01[1] * v02[0] - v01[0] * v02[1]); // left-handed coordinate (because pixel y-coordinate is going down)
}


/***
 * compute number of intersection of a ray against a line segment
 * @param org ray origin
 * @param dir ray direction (unit normal)
 * @param ps one of the two end points
 * @param pe the other end point
 * @return number of intersection
 */
int number_of_intersection_ray_against_edge(
    const Eigen::Vector2f &org,
    const Eigen::Vector2f &dir,
    const Eigen::Vector2f &ps,
    const Eigen::Vector2f &pe) {
  auto a = area(org, org + dir, ps);
  auto b = area(org, pe, org + dir);
  auto c = area(org, ps, pe);
  auto d = area(dir + ps, ps, pe);

  if (a * b > 0.f && d * c < 0.f) { return 1; }
  return 0;
  // the following code was a bug
  //auto d = area(org + dir, ps, pe);
  //if (a * b > 0.f && d * c > 0.f && fabs(d) > fabs(c)) { return 1; }
}


float func(float c2, float c1, float c0, float t) {
    return c2*t*t + c1*t + c0;
}


float dfunc(float c2, float c1, float t) {
    return 2*c2*t + c1;
}


// Actually, roots formula for quadratic function is easier than those two methods :)
float bisection(float c2, float c1, float c0, float left, float right) {
    float l = 0, r = 1;
    float mid;
    while (abs(left-right) > eps) {
        mid = (l + r) / 2;
        float tmp = func(c2, c1, c0, mid);
        if (tmp < eps) {
            break;
        }
        else if (tmp * func(c2, c1, c0, r) > 0) {
            r = mid;
        }
        else {
            l = mid;
        }
    }
    return mid;

}


float newton(float c2, float c1, float c0, float t0) {
    float t1;
    int itermax = 400, cnt = 0;
    while (cnt < itermax) {
        t1 = t0 - func(c2, c1, c0, t0) / dfunc(c2, c1, t0) / 2;
        cnt++;
        if (abs(t0-t1) < eps)
            break;
        else
            t0 = t1;
    }
    return t1;
}


int check(float c2, float c1, float c0,
          const Eigen::Vector2f &org,
          const Eigen::Vector2f &dir,
          const Eigen::Vector2f &ps,
          const Eigen::Vector2f &pc,
          const Eigen::Vector2f &pe) {
    int cnt = 0;
    float delta = sqrt(c1*c1-4*c2*c0);
    float t1 = (-c1+delta)/2/c2;
    if (t1 > 0 && t1 < 1) {
        float s1 = (ps[0]*(1-t1)*(1-t1) + pc[0]*2*(1-t1)*t1 + pe[0]*t1*t1 - org[0]) * dir[0];
        if (s1 > 0)
            cnt++;
    }
    float t2 = (-c1-delta)/2/c2;
    if (t2 > 0 && t2 < 1 && abs(t2-t1) > eps)  {
        float s2 = (ps[0]*(1-t2)*(1-t2) + pc[0]*2*(1-t2)*t2 + pe[0]*t2*t2 - org[0]) * dir[0];
        if (s2 > 0)
            cnt++;
    }
    return cnt;

}


/***
 *
 * @param org ray origin
 * @param dir ray direction (unit vector)
 * @param ps one of the two end points
 * @param pc control point
 * @param pe the other end point
 * @return the number of intersections
 */
int number_of_intersection_ray_against_quadratic_bezier(
    const Eigen::Vector2f &org,
    const Eigen::Vector2f &dir,
    const Eigen::Vector2f &ps,
    const Eigen::Vector2f &pc,
    const Eigen::Vector2f &pe) {
  // comment out below to do the assignment
  // return number_of_intersection_ray_against_edge(org, dir, ps, pe);
  // write some code below to find the intersection between ray and the quadratic
  float c2 = (ps[0] - 2*pc[0] + pe[0]) * (-dir[1]) + (ps[1] - 2*pc[1] + pe[1]) * dir[0];
  float c1 = (-2*ps[0] + 2*pc[0]) * (-dir[1]) + (-2*ps[1] + 2*pc[1]) * dir[0];
  float c0 = (ps[0] - org[0]) * (-dir[1]) + (ps[1] - org[1]) * dir[0];
  if (c2 < 0) {
      c2 = -c2;
      c1 = -c1;
      c0 = -c0;
  }
  float minv = c0 - c1*c1/4/c2;
  if (minv > 0) {
      return 0;
  }
  else {
      return check(c2, c1, c0, org, dir, ps, pc, pe);
  }

}

int main() {
  const auto input_file_path = std::filesystem::path(PROJECT_SOURCE_DIR) / ".." / "asset" / "r.svg";
  const auto [width, height, shape] = acg::svg_get_image_size_and_shape(input_file_path);
  if (width == 0) { // something went wrong in loading the function
    std::cout << "file open failure" << std::endl;
    abort();
  }
  const std::vector<std::string> outline_path = acg::svg_outline_path_from_shape(shape);
  const std::vector<std::vector<acg::Edge>> loops = acg::svg_loops_from_outline_path(outline_path);
  //
  std::vector<unsigned char> img_data(width * height, 255); // grayscale image initialized white
  for (unsigned int ih = 0; ih < height; ++ih) {
    for (unsigned int iw = 0; iw < width; ++iw) {
      const auto org = Eigen::Vector2f(iw + 0.5, ih + 0.5); // pixel center
      const auto dir = Eigen::Vector2f(60., 20.); // search direction
      int count_cross = 0;
      for (const auto &loop: loops) { // loop over loop (letter R have internal/external loops)
        for (const auto &edge: loop) { // loop over edge in the loop
          if (edge.is_bezier) { // in case the edge is a quadratic Bézier
            count_cross += number_of_intersection_ray_against_quadratic_bezier(
                org, dir,
                edge.ps, edge.pc, edge.pe);
          } else { // in case the edge is a line segment
            count_cross += number_of_intersection_ray_against_edge(
                org, dir,
                edge.ps, edge.pe);
          }
        }
      }
      if (count_cross % 2 == 1) { // Jordan's curve theory
        img_data[ih * width + iw] = 0; // paint black if it is inside
      }
    }
  }
  const auto output_file_path = std::filesystem::path(PROJECT_SOURCE_DIR) / "output.png";
  stbi_write_png(output_file_path.string().c_str(), width, height, 1, img_data.data(), width);
}
