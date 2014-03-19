#include <iostream>
#include <unistd.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>

#include "Imagine/Optim.h"
#include "parse_params.h"
#include "progressbar.h"
#include "stereo.h"

int main(int argc, char **argv)
{
  bool verbose = false;
  char *path = new char[100];
  strcpy(path,"./");  

  int nb_iter = 2;    
  int opt;   
  while((opt = getopt(argc, argv, "vp:n:")) > 0)
    switch(opt)
    {
      case 'v':
      verbose = true;
      break;
      case 'p':
      delete[] path;
      path = new char[strlen(optarg)+1];
      strcpy(path,optarg);
      break;    
      case 'n':
      nb_iter = strtol(optarg, NULL,10);
      break;      
      default:
      printf("La syntaxe est la suivante: ./PhotoStereo [-v] [-p path] [-n nb_iter]\n");
      return -1;
    }
  parse_ini(path);
  Camera::s = vec(LIGHT_X, LIGHT_Y, LIGHT_Z);
  if(nb_iter < 1) nb_iter = 1;
  if(nb_iter > 4) nb_iter = 4;  
  
  if(verbose)
    cout << "Verbose mode activated...\n";  
  
  camvec cl;
  featvec fl;
  cout << "Loading cameras...                    "; fflush(stdout);
  if(load_cam(verbose, path, cl, fl)) return 1; 
  cout << "done.\n";
  
  if (cl.size() < 4)
	{
		cerr << "\nThe PhotoStereo algorithm require at least 5 valid cameras... (only " << cl.size() << " valid cameras detected)\n";
		return 1;
  }

  cout << "Computing space search limits...      "; fflush(stdout);  
  SearchSpace ss(cl,fl);
  cout << "done.\n";
   
  if (ss.n_cam < 4)
	{
		cerr << "\nThe PhotoStereo algorithm require at least 5 valid cameras... (only " << ss.n_cam << " valid cameras detected)\n";
		return 1;
  }

  if(verbose)
  {
    cout << "number of cameras = " << ss.cl.size() << endl;  
    cout << "number of feature points = " << fl.size() << endl;    
  
    for(int i=0; i<ss.n_cam; i++)
      cout << "Cam " << i << "  Center: " << ss.cl[i]->center << "(norm: " << norm(ss.cl[i]->center) << ")    \tAxis: " << ss.cl[i]->axis << endl;  
  
    cout << "k: " << Camera::k << endl;
    cout << "I: " << ss.I << endl;
    cout << "J: " << ss.J << endl;
    cout << "K: " << ss.K << endl;  
    cout << "xmin: " << ss.xmin << endl;
    cout << "ymin: " << ss.ymin << endl;    
    cout << "zmin: " << ss.zmin << endl;
    cout << "xmax: " << ss.xmax << endl;
    cout << "ymax: " << ss.ymax << endl;        
    cout << "zmax: " << ss.zmax << endl;
    cout << "s: " << Camera::s << endl;    
  }
  
  if(false)
  {
    int width = 160;
    int height = 120;
    int n_j = 20;
    int cur_j = 10;  
    vec step = ss.get_step(width,height, n_j);
    Window mixwin = OpenWindow(width,height);
    Window Ewin = OpenWindow(width,height);    

    Event ev;
    camvec disp_cam;
    disp_cam = ss.cl;

    //disp_cam.push_back(ss.cl[0]); disp_cam.push_back(ss.cl[1]); disp_cam.push_back(ss.cl[2]); disp_cam.push_back(ss.cl[3]); 

    SetActiveWindow(mixwin);               
    display_plane(ss, width, height, ss.zmin+double(cur_j)*step[2], disp_cam);
    
    do {        
      cout << "j = " << cur_j << endl;  
      SetActiveWindow(mixwin);        
      display_plane(ss, width, height, ss.zmin+double(cur_j)*step[2], disp_cam);
      bool quit = false;
      do {
        GetEvent(500,ev); 
        switch (ev.type){
          case KEY:
          if(ev.key == 315 && cur_j<n_j-1) 
          {
            cur_j++;
            quit = true;      
          }
          if(ev.key == 317 && cur_j>0) 
          {
            cur_j--;        
            quit = true;
          }
          break;
          case CLICK:
          {
            Carac c;
            camvec cam;
            Camera *cam_ptr[4];
            vec p = (ss.xmin + double(ev.pix[0])*step[0]) * ss.I +
                    (ss.ymin + double(ev.pix[1])*step[1]) * ss.J +
                    (ss.zmin + double(cur_j)*step[2]) * ss.K;
            for(int i=0,k=0; i<ss.n_cam; i++)
            {
              vec P = ss.cl[i]->world2pict(p);
              if(!ss.cl[i]->is_out(P))
              {
                cam.push_back(ss.cl[i]);
                if(k<4)
                {
                  cam_ptr[k] = ss.cl[i];
                  k++;
                }
                cout << "color: " << ss.cl[i]->pix_grey(P) << endl;
              }
            }
            cout << "Point: " << p << endl;
            for(int i=0; i<4; i++)
            {
              cam_ptr[i]->precompute_info(p);
              //cout << "Pos cam: " << Camera::k*ss.world2plane(cam[i]->center-p) << endl;
              /*cout << "Point cam: " << ss.world2plane(cam_ptr[i]->world2cam(p)) << endl;
              vec l = Camera::s-cam_ptr[i]->world2cam(p);
              vec approx = ss.k*(cam[i]->center-p);
              approx /= pow(norm(approx),3.);
              cout << "CLV: " << ss.world2plane(l) << endl;
              cout << "approx LV: " << ss.world2plane(approx) << endl;
              cout << "LV: " << ss.world2plane(cam_ptr[i]->info.light_vector) << endl;            */
              //cout << "LV: " << cam_ptr[i]->info.light_vector << endl;                        
            }
            c = Carac(cam_ptr);
            cout << "n = " << c.n*ss.I << " " << c.n*ss.J << " " << c.n*ss.K << endl;
            cout << "p = " << c.p << endl;        
            cout << "a = " << c.a << endl;
            for(int i=0; i<4; i++)                            
              cout << "res = " << c.get_image_residual(cam_ptr[i]) << endl;
            double E = Ep(cam, p, c);
            cout << "E = " << E << endl;
            cout << "n = " << c.n*ss.I << " " << c.n*ss.J << " " << c.n*ss.K << endl;
            cout << "p = " << c.p << endl;        
            cout << "a = " << c.a << endl;        
          }
          break;
        }
      } while(!quit);
    } while (true);
  }
  
  cout << "Estimating depth and normals:\n";
  int W = 160;
  int H = 120;
  int D = 20; 
  int size_factor = 1;
  int depth_factor = 2;
  Image<pair<int,int> > depth_limit(W,H);  
  Image<int> depth;
  Image<double> smooth_depth;  
  Image<Carac> caracs;  
  for(int i=0; i<W; i++)
    for(int j=0; j<H; j++)  
    {
      depth_limit(i,j).first = 0;
      depth_limit(i,j).second = D-1;
    }
    
  for(int n = 0; n<nb_iter; n++)
  {  
    cout << "Iteration " << n+1 << " of " << nb_iter << ":\n"; fflush(stdout);
    graphcut_optimize(verbose, ss, D, depth_limit, depth, caracs);        
    
    cout << "Smoothing surface:\n"; fflush(stdout);  
    lsqr_optimize(verbose, ss, D, depth, caracs, smooth_depth);
    
    if(n<nb_iter-1)
    {
      D *= depth_factor;
      depth_limit = extrapolate(ss, smooth_depth, size_factor, depth_factor, D);     
    }     
  }
  
  if(true)
  {
    SetActiveWindow(OpenWindow(depth.width(),depth.height(),"Depth"));
    display(grey(depth));
    
    Image<Color> normals(depth.width(),depth.height());
    for(int u=0; u<depth.width(); u++)
      for(int v=0; v<depth.height(); v++)
      {
        vec n = ss.world2plane(caracs(u,v).n);
        double R = 255.*(n[0]+1.)/2.;
        double G = 255.*(n[1]+1.)/2.;
        double B = 255.*(n[2]+1.)/2.;                
        normals(u,v) = Color(byte(R), byte(G), byte(B));
      }
    SetActiveWindow(OpenWindow(normals.width(),normals.height(),"Normals"));      
    display(normals);      
    
    Image<Color> albedo(depth.width(),depth.height());
    for(int u=0; u<depth.width(); u++)
      for(int v=0; v<depth.height(); v++)
        albedo(u,v) = Color(byte(caracs(u,v).p[0]<0?0:(caracs(u,v).p[0]>255?255:caracs(u,v).p[0])), 
                            byte(caracs(u,v).p[1]<0?0:(caracs(u,v).p[1]>255?255:caracs(u,v).p[1])), 
                            byte(caracs(u,v).p[2]<0?0:(caracs(u,v).p[2]>255?255:caracs(u,v).p[2])));
    SetActiveWindow(OpenWindow(albedo.width(),albedo.height(),"Albedo"));      
    display(albedo);

    Image<Color> ambient(depth.width(),depth.height());
    for(int u=0; u<depth.width(); u++)
      for(int v=0; v<depth.height(); v++)
        ambient(u,v) = Color(byte(caracs(u,v).a[0]<0?0:(caracs(u,v).a[0]>255?255:caracs(u,v).a[0])), 
                             byte(caracs(u,v).a[1]<0?0:(caracs(u,v).a[1]>255?255:caracs(u,v).a[1])), 
                             byte(caracs(u,v).a[2]<0?0:(caracs(u,v).a[2]>255?255:caracs(u,v).a[2])));
    SetActiveWindow(OpenWindow(ambient.width(),ambient.height(),"Ambient"));      
    display(ambient);      

    SetActiveWindow(OpenWindow(smooth_depth.width(),
                               smooth_depth.height(),"smooth depth"));
    display(grey(smooth_depth));    
  }
  
  display_sphere();
  display_mesh(ss, smooth_depth);
}
/*----------------------------------------------------------------------------*/
Image<pair<int, int> > extrapolate(SearchSpace& ss, Image<double> depth, int size_factor, int depth_factor, int Dmax)
{
  const int W = depth.width();
  const int H = depth.height();
  vec step = ss.get_step(W,H,Dmax);
  double &Dz = step[2];
  Image<pair<int, int> > depth_limit(W*size_factor, H*size_factor);
  
  int size_nb = 30;
  for(int u=0; u<W; u++)
    for(int v=0; v<H; v++)
    {
      if(u==0 || v==0 || u==W-1 || v==H-1)
      {
        depth_limit(u,v).first = 0;
        depth_limit(u,v).second = Dmax-1;        
      }
      else
      {
        double zmin = depth(u,v);
        double zmax = depth(u,v);
        for(int x=u-size_nb; x<=u+size_nb; x++)
        {
          if(x<0 || x>=W) continue;
          for(int y=v-size_nb; y<=v+size_nb; y++)
          {
            if(y<0 || y>=H) continue;
            if(depth(x,y)<zmin) zmin = depth(x,y);
            if(depth(x,y)>zmax) zmax = depth(x,y);
          }
        }
        int jmin = (zmin-ss.zmin)/Dz;
        int jmax = (zmax-ss.zmin)/Dz + 0.99;
        int margin = 3;
        jmin = jmin-margin;
        jmax = jmax+margin;
        if(jmin<0) jmin = 0;
        if(jmax>=Dmax) jmax = Dmax-1;
        pair<int,int> limit(jmin,jmax);
        for(int x=size_factor*u; x<size_factor*W; x++)
          for(int y=size_factor*v; y<size_factor*H; y++)
            depth_limit(x,y) = limit;
      }
    }
  return depth_limit;
}
/*----------------------------------------------------------------------------*/
SearchSpace::SearchSpace(camvec& _cl, featvec fl)
{
  cl = _cl;
  n_cam = cl.size();
  
  vec center;
  ModifiedRansac<2>(cl.begin(), cl.end(), CenterEstimator(), CenterResidual(), center, tan(M_PI/5.)); 
 
  // I, J, K
  I = vec(0.,0.,0.);
  K = vec(0.,0.,0.);

  for(int i=0; i<n_cam; i++)
  {
    vec to_center = center-cl[i]->center;
    if(cl[i]->axis * to_center <= 0.) 
      continue; 
    else 
      to_center /= norm(to_center);
    K -= to_center;
    for(int j=i+1; j<n_cam; j++)
      I += cl[j]->center - cl[i]->center; 
  }
  
  if(norm(K) == 0.) K = cl.front()->axis;
  K.normalize();
  J = K ^ I;
  for(int i=0; i<n_cam && norm(J) == 0.; i++)
    J = K ^ cl[i]->center;
  J.normalize();
  I = J ^ K;
  
  
  for(int i=0; i<cl.size(); i++)
  {
    vec to_center = center-cl[i]->center;  
    if(cl[i]->axis * to_center <= 0. || to_center * K >= 0.)
    {
      cl.erase(cl.begin()+i);
      i--;
    }
  }  
  n_cam = cl.size();
  
  if(n_cam < 4)
    return;

  double z_step;
  double z0;
  const double filter_th = 0.2;
  const double enable_correl_th = 0.2;
  const double disable_correl_th = 0.2;  
  bool started = false;
  bool enabled = false;
  bool disabled = false;  
  double z;
  double c = 0.,last_c;  
   
  z = z0 = center*K;
  do
  {
    if(started)
      z += z_step;
    else
      z -= 1.;
      
    polygone inter;
    cl[0]->compute_clipping(I,J,K,z,inter);
    surface surf = inter.get_surface();
    double min_area = surf.area;
    for(int i=1; i<n_cam; i++)
    {
      polygone P;
      cl[i]->compute_clipping(I,J,K,z,P);
      inter = inter * P;
      if(P.get_surface().area < min_area)
        min_area = P.get_surface().area;        
    }
    last_c = c;
    c = inter.get_surface().area / min_area;
    
    if(started)
    {
      if(enabled && !disabled && c<=enable_correl_th) enabled = false;
      if(!enabled && c>enable_correl_th) {enabled = true; zmin = z;}
      if(enabled && c>disable_correl_th) disabled = true;
    }
    else if(c<filter_th || (last_c != 0. && c/last_c>0.9))
    {
      started = true;
      z_step = (z0-z) / 100;
    } 
    //cout << "z = " << z << " c = " << c << "    r = " << (last_c==0.?0.:c/last_c) << endl;
  }
  while((disabled && c>disable_correl_th) || 
        (!disabled && enabled && c>enable_correl_th) ||
        !enabled);
  zmax = z;   
  
  // xmin, xmax, ymin, ymax
  xmin = +std::numeric_limits<double>::infinity();
  xmax = -std::numeric_limits<double>::infinity();
  ymin = +std::numeric_limits<double>::infinity();
  ymax = -std::numeric_limits<double>::infinity();    
  z_step = (zmax - zmin) / 10;  
  for(int j=0; j<20; j++)
  {
    z = zmin + j*z_step;
    polygone P,inter;
    cl[0]->compute_clipping(I,J,K,z,inter);
    for(int i=1; i<n_cam; i++)
    {
      cl[i]->compute_clipping(I,J,K,z,P);
      inter = inter * P;
    }  
    if(!isinf(inter.ClipLT[0]) && inter.ClipLT[0] < xmin) xmin = inter.ClipLT[0];
    if(!isinf(inter.ClipLT[1]) && inter.ClipLT[1] < ymin) ymin = inter.ClipLT[1];
    if(!isinf(inter.ClipRB[0]) && inter.ClipRB[0] > xmax) xmax = inter.ClipRB[0];
    if(!isinf(inter.ClipRB[1]) && inter.ClipRB[1] > ymax) ymax = inter.ClipRB[1];     
  }
  
  // k
  Carac carac;
  double k;
  double k_inf = 0.;
  double k_sup = 5.;  
  
  double min_E;
  double min_k;
  
  double E;
  double step_k = 1.;
  
/*  do
  {
    min_E = std::numeric_limits<double>::infinity();
    for(k=k_inf; k<k_sup; k+=step_k)
    {
      if(k == 0) continue;
      Camera::k = k;
      E = 0;
      for(featiter feat = fl.begin(); feat != fl.end(); feat++)
      {
        E += Ep((*feat)->cam, (*feat)->pos, carac, 0.3);  
      }
      if(E < min_E)
      {
        min_E = E;
        min_k = k;
      }
    }
    k_inf = min_k - step_k;
    k_sup = min_k + step_k;
    if(k_inf<0.) k_inf = 0.;
    step_k /= 5.;
  }
  while(step_k>0.001);
  k = (min_k == 0.)?1.:min_k;
  Camera::k = k; */
}
/*----------------------------------------------------------------------------*/
vec SearchSpace::get_step(int W, int H, int D) const
{
  return vec((xmax-xmin)/W, (ymax-ymin)/H, (zmax-zmin)/D);
}
/*----------------------------------------------------------------------------*/














/*----------------------------------------------------------------------------*/
void display_plane(const SearchSpace& ss, int W, int H, double z, camvec& cam, bool draw_inside, Image<RGB<double> > *img)
{
  Image<RGB<double> > *iptr;
  if(img == NULL)
    iptr = new Image<RGB<double> >(W,H);  
  else
    iptr = img;
  vec step = ss.get_step(W,H,1);

  if(draw_inside)
    for(int u=0; u<W; u++)
      for(int v=0; v<H; v++)
      {
        (*iptr)(u,v) = RGB<double>(0.,0.,0.);
        for(int i=0; i<cam.size(); i++)
        {
          vec p = (ss.xmin+double(u)*step[0]) * ss.I +
                  (ss.ymin+double(v)*step[1]) * ss.J + z * ss.K;
          p = cam[i]->world2pict(p);
          if(!cam[i]->is_out(p))
            (*iptr)(u,v) += cam[i]->pix_color(p); 
        }
        (*iptr)(u,v) /= cam.size();
      }  
  display(color(*iptr));  
  
  polygone *P = new polygone[cam.size()];
  for(int i=0; i<cam.size(); i++)
  {
    cam[i]->compute_clipping(ss.I,ss.J,ss.K,z,P[i]); 
    P[i].draw(ss.xmin, ss.xmax, ss.ymin, ss.ymax, W, H, Color(255,0,0));
  }
  
  polygone inter = P[0];
  for(int i=1; i<cam.size(); i++)
    inter = inter*P[i];
  inter.draw(ss.xmin, ss.xmax, ss.ymin, ss.ymax, W, H, Color(0,255,0), 2);
 
  delete[] P;
  if(img == NULL) delete iptr;
}
/*----------------------------------------------------------------------------*/
void display_mesh(const SearchSpace &ss,const Image<double>& depth)
{
  double fact = 1.;
  int W = depth.width();
  int H = depth.height();
  Image<double> blured_depth(blur(depth,5.));  
  vec step = ss.get_step(W,H,1);  
	SetActiveWindow(OpenWindow(512,512,"3D",WindowVTK));
  int best_cam;
  double max_align = -1.;
  for(int i=0; i<ss.n_cam; i++)
    if(fabs(ss.K*ss.cl[i]->axis) > max_align)
    {
      best_cam = i;
      max_align = fabs(ss.K*ss.cl[i]->axis);
    }  

	Array<fPoint3> p(W*H);
	Array<Triangle> t(2*(W-1)*(H-1));
	Array<Color> col(W*H);
	if(false)   // use all views to paint the surface
	{
	  for (int j=0;j<H;j++)
		  for (int i=0;i<W;i++) { 
			  p[i+W*j]=fPoint3(double(i-W/2)*step[0],double(j-H/2)*step[1],blured_depth(i,j));
        vec P = (ss.xmin+double(i)*step[0]) * ss.I +
                (ss.ymin+double(j)*step[1]) * ss.J +
                blured_depth(i,j) * ss.K;
        RGB<double> c(0.,0.,0.);
        double n=0;        
        for(int k=0; k<ss.n_cam; k++)        
        {
          vec P2 = ss.cl[k]->world2pict(P);
          if(!ss.cl[k]->is_out(P2)) 
          {
            double dZdx = (blured_depth(i+1<W?i+1:i,j)-blured_depth(i-1<0?0:i-1,j))/2.;
            double dZdy = (blured_depth(i,j+1<H?j+1:j)-blured_depth(i,j-1<0?0:j-1))/2.;
            vec surf_n = ss.K - dZdx * ss.I - dZdy * ss.J;
            surf_n /= norm(surf_n);
            vec ray = ss.cl[k]->center-P;
            if(ray*surf_n > 0)
            {
              double coeff = ray*surf_n / norm(ray);
              n+=coeff;
              Color color = ss.cl[k]->pix_color(P2);
              c[0] += coeff*double(color[0]);
              c[1] += coeff*double(color[1]);
              c[2] += coeff*double(color[2]);
            }
          }
			  }
			  if(n != 0.)
			  {
			    c[0]/=n;
			    c[1]/=n;
			    c[2]/=n;			  			  
			  }
        col[i+W*j]=Color(byte(c[0]),byte(c[1]),byte(c[2]));			
		  }
  }
  else   // use a camera in front of the surface
  {
    for (int j=0;j<H;j++)
		  for (int i=0;i<W;i++) { 
			  p[i+W*j]=fPoint3(double(i-W/2)*step[0],double(j-H/2)*step[1],blured_depth(i,j));
        vec P = (ss.xmin+double(i)*step[0]) * ss.I +
                (ss.ymin+double(j)*step[1]) * ss.J +
                blured_depth(i,j) * ss.K; 	
        P = ss.cl[best_cam]->world2pict(P);
        if(ss.cl[best_cam]->is_out(P)) 
          col[i+W*j]=Color(0,0,0);
        else
        {  
          RGB<double> c = ss.cl[best_cam]->pix_color(P);
          col[i+W*j]=Color(byte(c[0]),byte(c[1]),byte(c[2]));
        }
      }			
  }
		
	for (int j=0;j<H-1;j++)
		for (int i=0;i<W-1;i++) {
			t[2*(i+j*(W-1))]=Triangle(i+W*j,i+1+W*j,i+W*(j+1));
			t[2*(i+j*(W-1))+1]=Triangle(i+1+W*j,i+1+W*(j+1),i+W*(j+1));
		}
	MeshVTK M(p.data(),W*H,t.data(),2*(W-1)*(H-1),PointColor);
	M.SetColors(col.data());
	MeshVTK Mn(p.data(),W*H,t.data(),2*(W-1)*(H-1),ConstantColor,GuessNormals);
	ShowMesh(M);
	bool text=true;
	Event ev;
	while (true) {
		GetEvent(5,ev);
		if (ev.type==CLICK && ev.button==1) {
			if (text) {
				HideMesh(M,false);
				ShowMesh(Mn,NoClipping,false);
			} else {
				HideMesh(Mn,false);
				ShowMesh(M,NoClipping,false);
			}
			text=!text;
		} else if (ev.type==CLICK && ev.button==3)
			break;
	}
	Terminate();
}
/*----------------------------------------------------------------------------*/
void display_sphere()
{
	SetActiveWindow(OpenWindow(512,512,"Normals Map"));	
	
	Image<Color> n_map(512,512);
	int R = 200;
	int d = 500;
  for(int x=-256;x<256;x++)  
  	for(int y=-256;y<256;y++)
  	{
  	  double a = x*x+y*y+d*d;
  	  double b = -2.*d*d;
  	  double c = d*d - R*R;
  	  double delta = b*b - 4.*a*c;
  	  if(delta >= 0.)
  	  {  	  
  	    double lambda = (-b - sqrt(delta)) / (2*a);
  	    vec dir(x,y,-d);
  	    vec o(0,0,d);
  	    o = o + lambda * dir;
  	    o /= double(R);
  	    o = (o+vec(1.,1.,1.))/2.;               
        o *= 255.;
  	    Color color;
  	    color[0] = byte(o[0]);
  	    color[1] = byte(o[1]);
  	    color[2] = byte(o[2]);
  	    n_map(256+x,256-y-1) = color;
  	  }  
  	  else
  	  n_map(256+x,256-y-1) = Color(255,255,255);
  	}
  display(n_map);	
}
/*----------------------------------------------------------------------------*/
