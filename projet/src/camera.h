#ifndef CAMERAH
#define CAMERAH

#include <Imagine/Images.h>
#include "polygone.h"

using namespace std;
using namespace Imagine;

typedef struct point_info {
  bool out;
  RGB<double> pix_color;
  double pix_grey;
  vec light_vector;  
} point_info; 

class Camera
{  
  public:
  int w2,h2;  
  int w,h;  
  double f,k1,k2;
	mat R,RT;
	vec t;
  static vec s;
  static double k;
	
	vec axis;
	vec center;
		
	Image<Color> img;
  point_info info;

  Camera();
  void load_cam(const char* path,const mat& _R,const vec& _t,const double& _f,const double& _k1,const double& _k2);
  
  void compute_clipping(const vec &I, const vec &J, const vec &K, double z, polygone &poly);
  void precompute_info(const vec& p_world);
  
  vec light_vector(const vec& p_world);
  inline vec world2cam(const vec& p_world) {return k*(R * p_world + t);}
  vec world2pict(const vec& p_world);   
  
  bool is_out(const vec& p_pict) const;
  
  RGB<double> pix_color(const vec& p_pict) const;
  double pix_grey(const vec& p_pict) const;
};

typedef std::vector<Camera*> camvec;
typedef camvec::iterator camiter;

typedef struct feature_point
{
  vec pos;
  Color color;  
  camvec cam;
  vector<vec> fp;  
} feature_point;

typedef std::vector<feature_point*> featvec;
typedef featvec::iterator featiter;

int load_cam(bool verbose, const string& dir, camvec& cl, featvec& fl);

#endif
