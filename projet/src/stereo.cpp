#include "stereo.h"
#include "progressbar.h"
#include "Imagine/Optim.h"
#include "maxflow/graph.h"
#include "lsqr.h"
#include "parse_params.h"

/*----------------------------------------------------------------------------*/
double Ep(camvec& cl, const vec& p_world, Carac& CaracEstimate, double OutRatio, double MinPb)
{
  camvec cam;  
  for(int i=0; i<cl.size();i++)
  {
    cl[i]->precompute_info(p_world);
    if(!cl[i]->info.out)
      cam.push_back(cl[i]);      
  }
  if(cam.size() < 4) return C0;
  
  int N = ModifiedRansac<4>(cam.begin(), cam.end(), CaracEstimator(), CaracResidual(), CaracEstimate, TAU, OutRatio, MinPb);
  if(N == 0) return C0; 

  double r, E = 0.;
  for(int i=0; i<cl.size(); i++)
  {
    r = CaracEstimate.get_image_residual(cl[i]);
    if(r < TAU)
      E += r;
  }

//  return (ETA * E - double(N)) / double(N);
//  return ETA * E / double(N) - double(N)/double(cl.size());
  return ETA * E / double(N) - double(N);
}
/*----------------------------------------------------------------------------*/
double En(const SearchSpace& ss, Image<vector<Carac> >& carac, Image<pair<int, int> > depth_limit, const vec step, Coords<3> M)
{
  vec& n = carac(M[0],M[1])[M[2]-depth_limit(M[0], M[1]).first].n; 
  if(norm(n) == 0.) return C1;

  vec p = (ss.xmin + double(M[0])*step[0]) * ss.I +
          (ss.ymin + double(M[1])*step[1]) * ss.J +
          (ss.zmin + double(M[2])*step[2]) * ss.K;
  double d = - n * p;
     
  double E = 0;
  int DzMax = Tj*sqrt(step[0]*step[0]+step[1]*step[1]) / step[2];
  
  list<pair<double, int> > dists;
  vector<vec> points;  
      
  for(int x = M[0]-1; x<=M[0]+1; x++)
  {
    if(x<0 || x>=depth_limit.width()) continue;
    for(int y = M[1]-1; y<=M[1]+1; y++)  
    {
      if(y<0 || y>=depth_limit.height()) continue;
      if(x == M[0] && y == M[1]) continue;

      int &first = depth_limit(x,y).first;
      int &last  = depth_limit(x,y).second;  

      dists.clear();
      points.resize(last-first+1);
      
      for(int j=M[2]-DzMax; j<=M[2]+DzMax; j++)
      {
        if(j<first || j>last) continue;
        if(norm(carac(x,y)[j-first].n) == 0.) continue;
        points[j-first] = 
             (ss.xmin + double(x)*step[0]) * ss.I +
             (ss.ymin + double(y)*step[1]) * ss.J +
             (ss.zmin + double(j)*step[2]) * ss.K;
        dists.push_back(pair<double, int>(points[j-first]*n + d,j)); 
      }
      
      if(dists.size() == 0)
        E += C1/8.;
      else
      {
        dists.sort();
        int j = dists.front().second;
        vec d = points[j-first] - p;
        E += (abs(j-M[2])+1) * fabs(n*d);
      }      
    }
  }
  return E;
}
/*----------------------------------------------------------------------------*/
void graphcut_optimize(bool verbose, SearchSpace& ss, int Dmax, Image<pair<int, int> >& depth_limit, Image<int>& depth, Image<Carac>& caracs)
{
  const int W = depth_limit.width();
  const int H = depth_limit.height();
  vec step = ss.get_step(W,H,Dmax);
  const double &Dx = step[0];
  const double &Dy = step[1];
  const double &Dz = step[2];
  double Emin = +std::numeric_limits<double>::infinity();
  double EpMean = 0;
  double EnMean = 0;
  double EsMean = 0;
  
  Image<vector<Carac> > carac(W, H);
  Image<vector<double> > E(W, H);
  depth = Image<int>(W,H);
  caracs = Image<Carac>(W, H);

  int maxD=0;
  double meanD=0;
  for(int u=0; u<W; u++)
  {
    int meanDH = 0;
    for(int v=0; v<H; v++)
    {
      E(u,v).resize(depth_limit(u,v).second-depth_limit(u,v).first+1);
      carac(u,v).resize(depth_limit(u,v).second-depth_limit(u,v).first+1);
          
      if(depth_limit(u,v).second-depth_limit(u,v).first+1>maxD)
        maxD = depth_limit(u,v).second-depth_limit(u,v).first+1;
      meanDH += depth_limit(u,v).second-depth_limit(u,v).first+1;
    }
    meanD += double(meanDH)/double(H);
  }
  meanD /= double(W);
  
  if(verbose)
    cout << "W = " << W << "   H = " << H << "   Dlim = " << Dmax << "   Dmax = " << maxD << "   Dmean = " << meanD << endl ;
   
  cout << "Computing Ep...                       ";
  ProgressBar PB_Ep(0,W); 
  for(int u=0; u<W; u++)
  {
    PB_Ep.update(u);
    for(int v=0; v<H; v++)
      for(int j=0; j<=depth_limit(u,v).second-depth_limit(u,v).first; j++)
      {
        vec p = (ss.xmin+double(u)*Dx) * ss.I + (ss.ymin+double(v)*Dy) * ss.J + (ss.zmin+double(j+depth_limit(u,v).first)*Dz) * ss.K;
        Coords<3> c(u,v,j);
        E(u,v)[j] = Ep(ss.cl, p, carac(u,v)[j], OutLierRatio, MinProba);         
      
        if(E(u,v)[j]<Emin)
          Emin = E(u,v)[j];
        EpMean += E(u,v)[j];
      }
  }
  EpMean /= W*H*meanD;
  EpMean -= Emin;
  PB_Ep.done();


  cout << "Computing En...                       ";    
  ProgressBar PB_En(0,W);   
  for(int u=0; u<W; u++)
  {
    PB_En.update(u); 
    for(int v=0; v<H; v++)
      for(int j=depth_limit(u,v).first; j<=depth_limit(u,v).second; j++)
      {
        double e = LAMBDA_N * En(ss, carac, depth_limit, step, Coords<3>(u,v,j));
        EnMean += e;
        E(u,v)[j-depth_limit(u,v).first] += e - Emin;
      }
  }
  EnMean /= W*H*meanD;
  PB_En.done();
  
  cout << "Bluring energy...                     ";  
  ProgressBar PB_blur(0,Dmax);   
  for(int j=0; j<Dmax; j++)
  {
    PB_blur.update(j);
    Image<double> E_blur(W,H);
    for(int x=0; x<W; x++)
      for(int y=0; y<H; y++)
        if(depth_limit(x,y).first<=j && j<=depth_limit(x,y).second)
        {
          E_blur(x,y) = E(x,y)[j-depth_limit(x,y).first];       
        }
        else
        {
          int n = 0;
          E_blur(x,y) = 0.;
          if(x>0 && y>0) {n++; E_blur(x,y)+=E_blur(x-1,y-1);}
          if(y>0) {n++; E_blur(x,y)+=E_blur(x,y-1);}          
          if(x>0) {n++; E_blur(x,y)+=E_blur(x-1,y);}                    
          if(n != 0)
            E_blur(x,y) /= n;  
        }
    E_blur = blur(E_blur,0.5);
    for(int x=0; x<W; x++)
      for(int y=0; y<H; y++)
        if(depth_limit(x,y).first<=j && j<=depth_limit(x,y).second)
          E(x,y)[j-depth_limit(x,y).first] = E_blur(x,y);
  }
  PB_blur.done();

  cout << "Computing graph cut...                ";   fflush(stdout);
  const int n = W*H*(Dmax+1);
  double infinity = std::numeric_limits<double>::infinity();
  Graph<double,double,double> G(n,2*n); 
  Image<int> IDS(W,H);

  for(int x=0; x<W; x++)
    for(int y=0; y<H; y++)
      IDS(x,y) = G.add_node(depth_limit(x,y).second - depth_limit(x,y).first + 2);
      
  for(int x=0; x<W; x++)
    for(int y=0; y<H; y++)
    {
      int D = depth_limit(x,y).second - depth_limit(x,y).first + 1;
      int id_first = IDS(x,y);
      int id_last  = IDS(x,y) + D; 
      
      G.add_tweights(id_first, infinity, 0.);
      G.add_tweights(id_last , 0., infinity);

      for(int j=depth_limit(x,y).first; j<=depth_limit(x,y).second; j++)
      {     
        double energy;
        if((x == 0 || y == 0 || x == W-1 || y == H-1) && j == depth_limit(x,y).first)
          energy = 0;
        else
          energy = E(x,y)[j-depth_limit(x,y).first];
        
        G.add_edge(IDS(x,y)+j-depth_limit(x,y).first,IDS(x,y)+j+1-depth_limit(x,y).first, energy, infinity);
      }
      
      if (x<W-1)
      { 
        int min = depth_limit(x,y).first<depth_limit(x+1,y).first?depth_limit(x,y).first:depth_limit(x+1,y).first;
        int max = depth_limit(x,y).second<depth_limit(x+1,y).second?depth_limit(x+1,y).second:depth_limit(x,y).second;        

        int id_first2 = IDS(x+1,y); 
        int id_last2  = IDS(x+1,y) + depth_limit(x+1,y).second - depth_limit(x+1,y).first + 1;         

        for(int j = min+1; j<max; j++)
        {
          int id1 = j<depth_limit(x,y).first?id_first:(
                    j>depth_limit(x,y).second?id_last:
                    IDS(x,y)+j-depth_limit(x,y).first);
          int id2 = j<depth_limit(x+1,y).first?id_first2:(
                    j>depth_limit(x+1,y).second?id_last2:
                    IDS(x+1,y)+j-depth_limit(x+1,y).first);
          double lambda;
          if(j == depth_limit(x,y).second+1 && j < depth_limit(x+1,y).first)
            lambda = LAMBDA_S*Dz*(depth_limit(x+1,y).first-j+1);          
          else if(j == depth_limit(x+1,y).second+1 && j < depth_limit(x,y).first)
            lambda = LAMBDA_S*Dz*(depth_limit(x,y).first-j+1);          
          else
            lambda = LAMBDA_S*Dz;
  			  G.add_edge(id1, id2, lambda, lambda);                 
        }
      }
      
		  if (y<H-1) 
      { 
        int min = depth_limit(x,y).first<depth_limit(x,y+1).first?depth_limit(x,y).first:depth_limit(x,y+1).first;
        int max = depth_limit(x,y).second<depth_limit(x,y+1).second?depth_limit(x,y+1).second:depth_limit(x,y).second;        

        int id_first2 = IDS(x,y+1);
        int id_last2  = IDS(x,y+1) + depth_limit(x,y+1).second - depth_limit(x,y+1).first + 1;

        for(int j = min+1; j<max; j++)
        {
          int id1 = j<depth_limit(x,y).first?id_first:(
                    j>depth_limit(x,y).second?id_last:
                    IDS(x,y)+j-depth_limit(x,y).first);
          int id2 = j<depth_limit(x,y+1).first?id_first2:(
                    j>depth_limit(x,y+1).second?id_last2:
                    IDS(x,y+1)+j-depth_limit(x,y+1).first);
          double lambda;
          if(j == depth_limit(x,y).second+1 && j < depth_limit(x,y+1).first)
            lambda = LAMBDA_S*Dz*(depth_limit(x,y+1).first-j+1);          
          else if(j == depth_limit(x,y+1).second+1 && j < depth_limit(x,y).first)
            lambda = LAMBDA_S*Dz*(depth_limit(x,y).first-j+1);          
          else
            lambda = LAMBDA_S*Dz;
  			  G.add_edge(id1, id2, lambda, lambda);                 
        }
      }		
  	}
  G.maxflow();
   
  for(int x=0; x<W; x++)
    for(int y=0; y<H; y++)
    {
      int id = IDS(x,y)+1;
			for(depth(x,y) = depth_limit(x,y).first; 
			    depth(x,y)<depth_limit(x,y).second && G.what_segment(id)==Graph<double,double,double>::SOURCE;
			    depth(x,y)++, id++) {}	    
		} 
  cout << "done." << endl;  

  if(verbose)
  {
    double EsMean = 0.;
    for(int x=0; x<W; x++)
      for(int y=0; y<H; y++)
      {
        if(x<W-1) EsMean += abs(depth(x,y)-depth(x+1,y))*LAMBDA_S*Dz;
        if(y<H-1) EsMean += abs(depth(x,y)-depth(x,y+1))*LAMBDA_S*Dz;        
      }
    EsMean /= W*H*meanD;
    cout << "Ep = " << EpMean << "   En = " << EnMean << "   Es = " << EsMean << endl;
  }
 
  for(int x=0; x<W; x++)
    for(int y=0; y<H; y++) 
      caracs(x,y) = carac(x,y)[depth(x,y)-depth_limit(x,y).first];
     
  //  show_energy_depth(ss, E, depth_limit, Dmax);   
/*
  bool done;
  int iter = 0;
  do
  {
    done = true;
    for(int x=0; x<W; x++)
      for(int y=0; y<H; y++) 
        if(norm(caracs(x,y).n) == 0.)
        {
          int n = 0;
          for(int u=x-1; u<=x+1; u++)
          {
            if(u<0 || u>=W) continue;
            for(int v=y-1; v<=y+1; v++)
            {
              if(v<0 || v>=H) continue;
              if(u==x && v==y) continue;
              if(norm(caracs(u,v).n) != 0.)
              {
                n++;
                caracs(x,y).n += caracs(u,v).n;
                caracs(x,y).a += caracs(u,v).a;
                caracs(x,y).p += caracs(u,v).p;
              }
            }
          }
          if(n!=0)
          {
            caracs(x,y).n /= norm(caracs(x,y).n);
            caracs(x,y).a /= n;
            caracs(x,y).p /= n;                    
          }
          else
            done = false;
        }
    iter++;
  }
  while(!done && iter<1);*/
}
/*----------------------------------------------------------------------------*/
void lsqr_optimize(bool verbose, SearchSpace& ss, int D, Image<int>& depth, Image<Carac>& caracs, Image<double>& optimal_depth)
{
  int W = depth.width();
  int H = depth.height();
  vec step = ss.get_step(W, H, D);
  double &Dx = step[0];
  double &Dy = step[1];
  double &Dz = step[2];    
  
  // Should be done by a real camera with depth interpolation?
  const double f = 1.;
  Image<long double> proj_depth(W,H);
  for(int x=0; x<W; x++)
    for(int y=0; y<H; y++)
      proj_depth(x,y) = double(depth(x,y))*Dz - (f + ss.zmax - ss.zmin);
  long double fx = double(W)*f / (ss.xmax-ss.xmin);
  long double fy = double(H)*f / (ss.ymax-ss.ymin);
  
  Surface surf(proj_depth, fx, fy, caracs); 
   
  cout << "Building system...                    ";
  System system = surf.allocate_system();
  cout << "done.\n";
  
  cout << "Solving system...                    ";
  Zvect optZ = solve_system(verbose, system, 0.001);
  cout << "done.\n";  
  
  optimal_depth = Image<double>(W,H);
  for(int x=0; x<W; x++)
    for(int y=0; y<H; y++)
      optimal_depth(x,y) = optZ[x*H+y] + f + ss.zmax;
}
/*----------------------------------------------------------------------------*/
void show_energy_depth(SearchSpace& ss, Image<vector<double> > E, Image<pair<int, int> >& depth_limit, int D)
{
  int W = E.width();
  int H = E.height();  
 	SetActiveWindow(OpenWindow(W,H, "Energy"));	
	vec step = ss.get_step(W,H,D);
	Image<double> energy(W,H);
	int j = 0;
	bool done = false;
	
	do
	{
    for(int x=0; x<W; x++)
      for(int y=0; y<H; y++)
        if(depth_limit(x,y).first<=j && j<=depth_limit(x,y).second)
          energy(x,y) = E(x,y)[j-depth_limit(x,y).first];
        else
          energy(x,y) = 0;
    display(grey(energy));
    bool quit = false;
    Event ev;
    do
    {
      GetEvent(500,ev); 
      switch (ev.type){
        case CLICK:
        if(ev.button != 1)
          done = quit = true;
        else
        {
          Carac c;
          camvec cam;
          Camera *cam_ptr[4];
          vec p = (ss.xmin + double(ev.pix[0])*step[0]) * ss.I +
                  (ss.ymin + double(ev.pix[1])*step[1]) * ss.J +
                  (ss.zmin + double(j)*step[2]) * ss.K;
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
          double energy;
          if(depth_limit(ev.pix[0],ev.pix[1]).first<=j && j<=depth_limit(ev.pix[0],ev.pix[1]).second)
            energy = E(ev.pix[0],ev.pix[1])[j-depth_limit(ev.pix[0],ev.pix[1]).first];
          else
            energy = 0.;
          cout << "E current = " << energy << endl;      
          cout << "jmin = " << depth_limit(ev.pix[0],ev.pix[1]).first << endl;
          cout << "jmax = " << depth_limit(ev.pix[0],ev.pix[1]).second << endl;          
          for(int j2=depth_limit(ev.pix[0],ev.pix[1]).first; j2<=depth_limit(ev.pix[0],ev.pix[1]).second; j2++)
            cout << "E(" << j2 << ") = " << E(ev.pix[0],ev.pix[1])[j2-depth_limit(ev.pix[0],ev.pix[1]).first] << endl;
        }
        break;
        case KEY:
        if(ev.key == 315 && j<D-1) 
        {
          j++;
          cout << "j = " << j << endl;   
          quit = true;     
        }
        if(ev.key == 317 && j>0) 
        {
          j--;        
          cout << "j = " << j << endl;
          quit = true;
        }
        break;    
      }
    }    
    while(!quit);  
	}
	while(!done);
}
