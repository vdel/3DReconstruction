#include <iostream>
#include <fstream>

#include "camera.h"

vec Camera::s = vec(0.,0.,0.);   
double Camera::k = 1.;

/*----------------------------------------------------------------------------*/
Camera::Camera()
{
}
/*----------------------------------------------------------------------------*/
void Camera::load_cam(const char* path,const mat& _R,const vec& _t,const double& _f,const double& _k1,const double& _k2)
{
  R  = _R;
  RT = transpose(R);
  t  = _t;
  f  = fabs(_f);
  k1 = _k1;
  k2 = _k2;
  load(img,path);
  h = img.height();
  w = img.width();
  h2 = h/2;
  w2 = w/2;
  center = -transpose(R)*t;
  axis = -vec(R.getrow(2).data());
}
/*----------------------------------------------------------------------------*/
vec Camera::light_vector(const vec& p_world)
{
  vec l = s - world2cam(p_world);
  return RT*l/pow(norm(l),3.);
}
/*----------------------------------------------------------------------------*/
vec Camera::world2pict(const vec& p_world)
{
  vec p = k*(R * p_world + t);           // camera coordinate
  p = -p/p[2];                     // perspective division
  p *= f;                          // pixel coordinates
//    double n2 = pow(norm(p), 2.);
//    p *= (1. + n2 * (k1 + n2 * k2));  // distortion correction
  p[0] = p[0] + w2;
  p[1] = h2 - p[1];
  return p;
}
/*----------------------------------------------------------------------------*/
void Camera::compute_clipping(const vec &I, const vec &J, const vec &K, double z, polygone &poly)
{
	vec dir[4] = { vec(w2, h2, -f), vec(-w2, h2, -f), vec(-w2, -h2, -f), vec(w2, -h2, -f) };
  
  double centerK = (z*K-center)*K;
  double axisK = axis*K;
    
	for(int i=0; i<4; i++)
    dir[i] = transpose(R) * dir[i];
    
  poly.begin_def();  
	for(int i=0; i<4; i++)
	{
	  double lambda;
    double dsk = dir[i]*K;
    if(dsk == 0.) lambda = -1.; else lambda = centerK/dsk;
    if(dsk*axisK<0.) // Towards infinity
    {
      vec n[2];
      n[0] = dir[(i+3)%4] ^ dir[i];
      n[1] = dir[i] ^ dir[(i+1)%4];
      
      for(int k=0; k<2; k++)
      {
        vec d = K ^ n[i];
        if(norm(d) != 0)
     	    poly.add_point(vec(d*I, d*J, 0.)); 
      }       	      
    }
    else
    {
      vec M = center+lambda*dir[i];	      
	    poly.add_point(vec(M*I, M*J, 1.));          
    }       
  }
  poly.end_def();    
}
/*----------------------------------------------------------------------------*/
void Camera::precompute_info(const vec& p_world)
{
  vec p_pict = world2pict(p_world);
  if(!is_out(p_pict))
  {
    info.out = false;
    info.pix_color = pix_color(p_pict);
    info.pix_grey = pix_grey(p_pict);
    info.light_vector = light_vector(p_world);
  }
  else
    info.out = true;
}
/*----------------------------------------------------------------------------*/
bool Camera::is_out(const vec& p_pict) const
{
  return p_pict[0]<0 || p_pict[1]<0 || p_pict[0]>=w || p_pict[1]>=h;
}
/*----------------------------------------------------------------------------*/
RGB<double> Camera::pix_color(const vec& p_pict) const
{
  return img.interpolate(p_pict[0], p_pict[1]);
}
/*----------------------------------------------------------------------------*/
double Camera::pix_grey(const vec& p_pict) const
{
  Color c = img.interpolate(p_pict[0], p_pict[1]);
  return 0.299*double(c[0]) + 0.587*double(c[1]) + 0.114*double(c[2]);
}
/*----------------------------------------------------------------------------*/









/*----------------------------------------------------------------------------*/
int load_cam(bool verbose, const string& dir, camvec& cl, featvec& fl)
{
  char txt[200];
  ifstream bundle;
  ifstream list;
  string img;
  string list_file = dir;
  list_file += "list.txt";
  string bundle_file = dir;
  bundle_file += "bundle/bundle.out";
  
  bundle.open(bundle_file.c_str());
	if (!bundle.is_open())
	{
		cerr << "\nFile " << bundle_file << " not found." << endl;
    cerr << "La syntaxe est la suivante: ./PhotoStereo [-v] [-p path] [-n nb_iter]\n";		
		return 1;
  }  bundle >> txt >> txt >> txt >> txt;
  
  int n_cam, n_points;
  bundle >> n_cam >> n_points;
  
  // Reading camera parameters
  mat *R  = new mat[n_cam];
  vec *t  = new vec[n_cam];
  double *f  = new double[n_cam];
  double *k1 = new double[n_cam];  
  double *k2 = new double[n_cam];    
  bool *is_valid = new bool[n_cam];
  
  for(int n=0; n<n_cam; n++)
  {
    bundle >> f[n] >> k1[n] >> k2[n];
    
    for (int r=0;r<3;r++)
      bundle >> R[n](r,0) >> R[n](r,1) >> R[n](r,2);
    bundle >> t[n][0] >> t[n][1] >> t[n][2];
    
    is_valid[n] = (t[n][0]!=0. || t[n][1]!=0. || t[n][2]!=0.);
  }
  
  // Reading pictures files
  list.open(list_file.c_str());
	if (!list.is_open())
	{
		cerr << "\nFile " << list_file << " not found." << endl;
		return 1;
  }
  
  vector<int> id_c;
  id_c.resize(n_cam);
  if(verbose)  
    cout << endl;
  for(int n=0; n<n_cam; n++)
  {
    list >> txt;
    if(is_valid[n])
    {
      img = dir;
      img += txt;
      id_c[n] = cl.size();                  
      cl.push_back(new Camera);
      cl.back()->load_cam(img.c_str(), R[n], t[n], f[n], k1[n], k2[n]);
      if(verbose)
        cout << "Loaded image '" << img.c_str() << "'.\n";
    }  
    else
      id_c[n] = -1;
    list >> txt;
    list >> txt; 
  }  
  list.close();
  
  // Reading feature points
  for(int n=0; n<n_points; n++)  
  {
    int num_cams;
    int R,G,B;
    fl.push_back(new feature_point);
    bundle >> fl.back()->pos[0]
           >> fl.back()->pos[1]
           >> fl.back()->pos[2];
    bundle >> R
           >> G
           >> B;
    fl.back()->color = Color(R,G,B);       
    bundle >> num_cams;           
    for(int m=0; m<num_cams; m++)
    {
      int i;
      vec p;
      bundle >> i >> txt >> p[0] >> p[1];
      
      if(id_c[i] == -1)
        continue;     
      
      p[0] += cl[id_c[i]]->w2; 
      p[1]  = cl[id_c[i]]->h2 - p[1];          
      p[2] = 1;      
      fl.back()->fp.push_back(p);
      fl.back()->cam.push_back(cl[id_c[i]]);      
    }
    if(fl.back()->cam.size()<5)
    {
      delete fl.back();
      fl.pop_back();
    }    
  }
  bundle.close();
  
  delete[] R;
  delete[] t;
  delete[] f;
  delete[] k1;
  delete[] k2;
  delete[] is_valid; 
  return 0;
}
/*----------------------------------------------------------------------------*/
