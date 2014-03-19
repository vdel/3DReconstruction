#ifndef POLYGONEH
#define POLYGONEH

#include <Imagine/Images.h>
#include <vector>

using namespace std;
using namespace Imagine;

typedef FVector<double,3> vec;
typedef FMatrix<double,3,3> mat;

class surface
{
  public:
  double area;
  double plane_part;
  
  surface() {area = 0.; plane_part = 0.;}
  surface(double _a, double _pp) {area = _a; plane_part = _pp;}
  surface operator + (const surface& s) {return surface(area+s.area, plane_part+s.plane_part);}
  surface operator / (const double& c) {return surface(area/c, plane_part/c);}  
  double operator / (const surface& s) 
  {
    if(s.plane_part != 0.) return plane_part/s.plane_part;
    else if(s.area != 0.)  return area/s.area;
                           return 0.;
  }
};

class segment
{
  public:
  vec A, B, n;
  segment(const segment& s);
  segment(vec& p1, vec& p2); 
  
  enum {LEFT, INTER, RIGHT};
  int intersect(const vec& n, int& side);  
};

class polygone
{
  protected:
  void clear();
  void clear_points();  
  void clear_edges();  
  void inter_halfspace(const vec& n);
     
  public:
  bool is_infinite;  
  vec ClipLT, ClipRB;
  vector<vec*> points;
  vector<segment*> edges;
      
  polygone();
  polygone(const polygone& p);
  void operator = (const polygone& p);
  
  void begin_def();
  void add_point(const vec& p);
  void end_def();
   
  polygone operator * (const polygone& P2);  // Intersection
  bool is_in(const vec& p);
 
  surface get_surface();
  
  void draw(double xmin, double xmax, double ymin, double ymax, int W, int H, Color c = Color(255,0,0), int width = 1);
};

#endif
