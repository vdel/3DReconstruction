#ifndef CARACH
#define CARACH

#include "camera.h"

class Carac
{
  public:
  vec n;  //normal
  vec p;  // albedo
  vec a;  // ambient
  
  Carac();
  Carac(Camera **cam);

  double get_image_residual(Camera* cam) const;  
};

/*----------------------------------------------------------------------------*/
// Functor: estimates normal, albedo and ambient from 4 pairs
class CaracEstimator {
public:	
	template <class tcaracs>
	void operator() (FArray<camiter,4>& Icam, tcaracs caracs) const {
	  Camera *cam[4];
	  for(int i=0; i<4; i++)
	    cam[i] = *(Icam[i]);
	  
	  Carac c(cam);
		if(norm(c.n) != 0.)
    	*caracs++ = c;
	}
};
/*----------------------------------------------------------------------------*/
// Functor: estimates the energy gi(p,j)
class CaracResidual {
public:	
	inline double operator() (const Carac& carac, Camera* view) const {
  	return carac.get_image_residual(view);
	}
};

#endif
