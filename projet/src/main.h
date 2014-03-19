#ifndef MAINH
#define MAINH

#include "camera.h"

/*----------------------------------------------------------------------------*/
// Functor: estimates center of the view given two camera
class CenterEstimator {
public:	
	template <class vectors>
	void operator() (FArray<camiter,2>& Icam, vectors vects) const {
	  Camera *cam[2];
	  for(int i=0; i<2; i++)
	    cam[i] = *(Icam[i]);
	  
	  vec k = cam[0]->axis ^ cam[1]->axis;
	  if(norm(k) != 0.)
	  {
	    k /= norm(k);
	    vec c0c1 = cam[1]->center - cam[0]->center;
  	  vec j = k ^ cam[1]->axis;
  	  double lambda = (c0c1 * j) / (cam[0]->axis * j);  	
  	  if(lambda>0)
  	  {  
  	    vec O = cam[0]->center + lambda * cam[0]->axis;    	    	    	  
  	    double dist = c0c1 * k;
  	    if((O + dist*k - cam[1]->center) * cam[1]->axis > 0)
  	      *vects++ = O + (dist/2.)*k;
  	  }
	  }
	}
};
/*----------------------------------------------------------------------------*/
// Functor: estimates the distance between the center and the axis
class CenterResidual {
public:	
	double operator() (const vec& center, Camera* view) const {
	  double lambda = (center - view->center) * view->axis;
	  if(lambda>0)
	  {
      vec p = center - (view->center + lambda*view->axis);
      return norm(p)/lambda;	  
    }
    else
      return std::numeric_limits<double>::infinity();
	}
};





/*----------------------------------------------------------------------------*/
class SearchSpace
{
  public:
  camvec cl;
  int n_cam;
  
  vec    I,  J,  K;   // base of planes for plane-sweep stereo
  double xmin, ymin, zmin;  // beginning of the search space
  double xmax, ymax, zmax;  // end of the search space
  
  SearchSpace(camvec& _cl, featvec fl);
  vec get_step(int W, int H, int D) const;
  inline vec world2plane(const vec& p) {return vec(p*I,p*J,p*K);}
  inline vec plane2world(const vec& p) {return p[0]*I+p[1]*J+p[2]*K;}  
};


/*----------------------------------------------------------------------------*/
Image<pair<int, int> > extrapolate(SearchSpace& ss, Image<double> depth, int size_factor, int depth_factor, int Dmax);

void display_plane(const SearchSpace& ss, int W, int H, double z, camvec& cam, bool draw_inside=true, Image<RGB<double> > *img=NULL);
void display_mesh(const SearchSpace &ss,const Image<double>& depth);
void display_sphere();

#endif
