#include "polygone.h"


/*----------------------------------------------------------------------------*/
segment::segment(const segment& s)
{
  n = s.n;
  A = s.A;
  B = s.B;
}
/*----------------------------------------------------------------------------*/
segment::segment(vec& p1, vec& p2)
{
  if(p1[2] == 0.)
  { 
    A = -p1;
    B = p2;
    n = p1;
    n = vec(-n[1], n[0], 0.); n[2] = -B*n;
  }
  else if(p2[2] == 0.)
  {
    A = p1;
    B = p2;
    n = p2;
    n = vec(-n[1], n[0], 0.); n[2] = -A*n;
  }
  else
  {
    A = p1;
    B = p2;
    n = B-A;
    n = vec(-n[1], n[0], 0.); n[2] = -A*n;   
  } 
}
/*----------------------------------------------------------------------------*/
int segment::intersect(const vec& N, int& side)
{
  double a = (A * N);
  double b = (B * N);
  double s = a * b;
  
  if(s < 0. || (s == 0. && a+b > 0.)) // intersection - updates the segment
  {
    if(a <= 0.)
    {
      A = N ^ n;
      A = A / A[2];
      side = LEFT;
    }  
    else
    {
      B = N ^ n;
      B = B / B[2];
      side = RIGHT;
    }      
    return INTER;
  }
  else        // no intersection - points on the same side
  {
    if(a+b > 0.) side = LEFT; else side = RIGHT;
    return side;
  }
}
/*----------------------------------------------------------------------------*/







/*----------------------------------------------------------------------------*/
polygone::polygone()
{
}
/*----------------------------------------------------------------------------*/
polygone::polygone(const polygone& p)
{
  points.resize(p.points.size());
  for(int i=0; i<p.points.size(); i++)
    points[i] = new vec(*p.points[i]);
    
  edges.resize(p.edges.size());  
  for(int i=0; i<p.edges.size(); i++)
    edges[i] = new segment(*p.edges[i]);

  is_infinite = p.is_infinite;
  ClipLT = p.ClipLT;    
  ClipRB = p.ClipRB;  
}
/*----------------------------------------------------------------------------*/
void polygone::operator = (const polygone& p)
{
  clear();
  points.resize(p.points.size());
  for(int i=0; i<p.points.size(); i++)
    points[i] = new vec(*p.points[i]);
    
  edges.resize(p.edges.size());  
  for(int i=0; i<p.edges.size(); i++)
    edges[i] = new segment(*p.edges[i]);
    
  is_infinite = p.is_infinite;
  ClipLT = p.ClipLT;    
  ClipRB = p.ClipRB;    
}
/*----------------------------------------------------------------------------*/     
void polygone::clear()
{
  clear_points();
  clear_edges();
}
/*----------------------------------------------------------------------------*/     
void polygone::clear_points()
{
  is_infinite = false;
  for(int i=0; i<points.size(); i++)
    delete points[i];
  points.clear();

  ClipLT[0] = +std::numeric_limits<double>::infinity();
  ClipLT[1] = +std::numeric_limits<double>::infinity();
  ClipRB[0] = -std::numeric_limits<double>::infinity();
  ClipRB[1] = -std::numeric_limits<double>::infinity();   
}
/*----------------------------------------------------------------------------*/
void polygone::clear_edges()
{  
  for(int i=0; i<edges.size(); i++)
    delete edges[i];  
  edges.clear();  
}
/*----------------------------------------------------------------------------*/
void polygone::begin_def()
{
  clear();
}
/*----------------------------------------------------------------------------*/
void polygone::end_def()
{
  const int n = points.size();
  for(int i=0; i<n; i++)
    if((*points[i])[2] != 0. || (*points[(i+1)%n])[2] !=0.)
    {
      edges.resize(edges.size()+1);
      edges[edges.size()-1] = new segment(*points[i], *points[(i+1)%n]);
    }  
}
/*----------------------------------------------------------------------------*/
void polygone::add_point(const vec& p)
{
  vec* P = new vec(p);
  if((*P)[2] == 0.) 
  {
    is_infinite = true;
    points.push_back(P);
  }
  else
  {
    *P = *P / (*P)[2];
    if((*P)[0] < ClipLT[0]) ClipLT[0] = (*P)[0];
    if((*P)[1] < ClipLT[1]) ClipLT[1] = (*P)[1];
    if((*P)[0] > ClipRB[0]) ClipRB[0] = (*P)[0];
    if((*P)[1] > ClipRB[1]) ClipRB[1] = (*P)[1];        
    points.push_back(P);
  }
}
/*----------------------------------------------------------------------------*/
polygone polygone::operator * (const polygone& P2)
{
  if(edges.size() == 0) return *this;
  else
  {
    polygone I = P2;
    for(int i=0; i<edges.size(); i++)
      I.inter_halfspace(edges[i]->n);
    return I;
  }
}
/*----------------------------------------------------------------------------*/
bool polygone::is_in(const vec& p)
{
  int i;
  for(i=0; i<edges.size(); i++)
    if(edges[i]->n*p<=0) break;
  return i == edges.size();
}
/*----------------------------------------------------------------------------*/
void polygone::inter_halfspace(const vec& n)
{
  int i,j = 0;
  int num_inter = 0;
  int side;
  int init_side;
  vec A, B;

  for(i=0; i<edges.size(); i++)
  {  
    switch(edges[i]->intersect(n,side))
    {
      case segment::LEFT: //segment is in the half space we want: we keep it
      edges[j] = edges[i];
      j++;
      break;
      case segment::INTER: //segment is cut: we keep its left part
      if(num_inter)
      {
        num_inter++;
        
        if(side == segment::LEFT)        
          B = edges[i]->A;
        else
        {
          B = edges[i]->B;        
          edges[j] = edges[i];
          j++;            
        }
          
        if((A[2] != 0. || B[2] !=0.) && A != B)
        {
          if(i <= j)
          {
            edges.resize(edges.size()+1);
            for(int k=edges.size()-1; i<k; k--)
              edges[k] = edges[k-1];
            i++;  
          }
          
          if(side == segment::LEFT)
          {
            edges[j] = new segment(A, B);
            j++;
            edges[j] = edges[i];
          }
          else
            edges[j] = new segment(B, A);
          j++;
        }
      }
      else
      {
        num_inter++;
        init_side = side;
        
        if(side == segment::LEFT)
          A = edges[i]->A;
        else 
          A = edges[i]->B;

        edges[j] = edges[i];
        j++;            
      }
      break;
      case segment::RIGHT: // segment is in the wrong half space: we delete it
      delete edges[i];
      break;
    }
  }
  if(num_inter == 1) // Current segment is an infinite edge
  {
    B = vec(n[1],-n[0],0.);
    if(B[2] != 0. || A[2] !=0.)
    {
      if(init_side == segment::LEFT)
        edges[j] = new segment(B, A);      
      else
        edges[j] = new segment(A, B);
      j++;
    }
  }  
  edges.resize(j);
  
  // recomputes points
  clear_points();    
  for(i=0; i<edges.size(); i++)
  {
    if(edges[i]->A[2] == 0.)
      add_point(-edges[i]->A);
    add_point(edges[i]->B);          
  }  
}
/*----------------------------------------------------------------------------*/
surface polygone::get_surface()
{
  if(is_infinite)
  {
    vec *I[2];
    int j = 0;
    for(int i=0; i<points.size(); i++)
    {
      if((*points[i])[2] == 0.)
      {
        I[j] = points[i];
        j++;
      }  
      if(j == 2) break;
    }  
    if(j != 2) return surface(0., 0.);
    else
    {
      double c = (*I[0])*(*I[1])/norm(*I[0])/norm(*I[1]);
      return surface(std::numeric_limits<double>::infinity(), acos(c)/M_PI/2.);      
    }  
  }
  else
  {
    double s = 0.;
    
    for(int i=0; i<(signed int)(points.size())-2; i++)
    {
      vec c1 = (*points[i+1])-(*points[i]);
      vec c2 = (*points[points.size()-1])-(*points[i]);
      s += fabs(c1[0]*c2[1]-c1[1]*c2[0])/2.;  
    }
    
    return surface(s, 0.);  
  }
}
/*----------------------------------------------------------------------------*/
void polygone::draw(double xmin, double xmax, double ymin, double ymax, int W, int H, Color c, int width)
{
  double Dx = (xmax-xmin) / W;
  double Dy = (ymax-ymin) / H;  
  if(points.size() == 0) return;
  int offset = 0;
  vector<int> x;
  vector<int> y;
  bool last_was_inf = (*points[points.size()-1])[2];
  for(int j=0; j<points.size(); j++)
  {
    if((*points[j])[2] == 0.)
    {
      vec p, dir;
      if(last_was_inf)
      {
        p = *points[(j+1)%points.size()];
        dir = vec((*points[j])[1],-(*points[j])[0],0.);
      }
      else
      {
        p = *points[(j+points.size()-1)%points.size()];
        dir = vec(-(*points[j])[1],(*points[j])[0],0.);
      }
      dir[2] = -p*dir;
      vec h,v;
      if((*points[j])[0]<=0) h = vec(1.,0.,-xmin)^dir;
                        else h = vec(1.,0.,-xmin*W*Dx)^dir;
      if((*points[j])[1]<=0) v = vec(0.,1.,-ymin)^dir;
                        else v = vec(0.,1.,-ymin*H*Dy)^dir;
      if(h[2] != 0.) h /= h[2];
      else {h[0] = -1; h[1] = -1;}
      if(v[2] != 0.) v /= v[2];
      else {v[0] = -1; v[1] = -1;}
      if(h[1]<0 || h[1]>=H)
      {
        x.push_back(v[0]);
        y.push_back(v[1]);
      }
      else
        x.push_back(h[0]);
        y.push_back(h[1]);      
      if(last_was_inf) 
        offset = x.size()-1;
      last_was_inf = true;
    }
    else
    {
      x.push_back(((*points[j])[0]-xmin)/Dx);
      y.push_back(((*points[j])[1]-ymin)/Dy);      
    }
  }
  int *cx = new int[points.size()];
  int *cy = new int[points.size()];
  for(int j=0; j<points.size(); j++)
  {
    cx[(j+points.size()-offset)%points.size()] = x[j];
    cy[(j+points.size()-offset)%points.size()] = y[j];
  }
  if(is_infinite)
    for(int j=1; j<points.size(); j++)
      DrawLine(cx[j-1],cy[j-1],cx[j],cy[j],c,width);
  else
    DrawPoly(cx,cy,points.size(),c,width);
  delete[] cx;
  delete[] cy;
}
/*----------------------------------------------------------------------------*/
