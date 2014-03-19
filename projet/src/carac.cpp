#include <iostream>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>

#include "carac.h"
#include "parse_params.h"

using namespace Imagine;

/*----------------------------------------------------------------------------*/
Carac::Carac()
{
  n = vec(0.,0.,0.);
  p = vec(0.,0.,0.);
  a = vec(0.,0.,0.);    
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Compute caracteritics from 4 views
/*----------------------------------------------------------------------------*/
// Normal version
Carac::Carac(Camera **cam)
{
  Matrix<double> A(4,4);
	Vector<double> B(4);
     
  for(int i=0; i<4; i++)
  {
    A(i,0) = cam[i]->info.light_vector[0];
    A(i,1) = cam[i]->info.light_vector[1];
    A(i,2) = cam[i]->info.light_vector[2];
    A(i,3) = 1.;
    B[i] = cam[i]->info.pix_grey/LIGHT_SOURCE;
  }

  A = inverse(A);
        
  if(norm(A) == 0.)  // out of a picture or A not invertible
  {
    for(int i = 0; i<3; i++)
    {
      n[i] = 0.;
      a[i] = 0.;
      p[i] = 0.;
    }
  }
  else
  {
    B = A*B;
    n[0] = B[0];
    n[1] = B[1];
    n[2] = B[2];
    n /= norm(n);
    const double a_grey = B[3]*LIGHT_SOURCE;
    p = vec(0.,0.,0.);
    
    for(int i=0; i<4; i++)
    {
      double r = cam[i]->info.light_vector*n+a_grey/255.;
      for(int c=0; c<3; c++)
        p[c] += cam[i]->info.pix_color[c] / r;
    }
    p /= 4.;
    for(int c=0; c<3; c++)
      a[c] = p[c] * a_grey / 255.;  
  }
}
/*----------------------------------------------------------------------------*/
double Carac::get_image_residual(Camera* cam) const
{
  double coeff = cam->info.light_vector*n;
  double r = 0.;  
  for(int c=0; c<3; c++)
    r += fabs(cam->info.pix_color[c] - p[c] * coeff - a[c]);
  return r;
}
/*----------------------------------------------------------------------------*/




















/* Compute caracteritics from 4 views
/*----------------------------------------------------------------------------*/
// Grey scale carac
/*Carac::Carac(Camera **cam, bool _v)
{ 
  verbose = _v;
  Matrix<double> A(4,4);
	Vector<double> B(4);
       
  for(int i=0; i<4; i++)
  {
    B[i]   = cam[i]->info.pix_grey/LIGHT_SOURCE;
    A(i,0) = cam[i]->info.light_vector[0];
    A(i,1) = cam[i]->info.light_vector[1];
    A(i,2) = cam[i]->info.light_vector[2];
    A(i,3) = 1.;
  }

  A = inverse(A);
        
  if(norm(A) == 0.)  // out of a picture or A not invertible
  {
    for(int i = 0; i<3; i++)
    {
      n[i] = 0.;
      a[i] = 0.;
      p[i] = 0.;
    }
  }
  else
  {
    B = A*B;
    n[0] = B[0];
    n[1] = B[1];
    n[2] = B[2];
    if(verbose)
      cout << B << endl;
    p[0] = norm(n);
    n /= p[0];
    a[0] = B[3]*LIGHT_SOURCE;
  }
}
//----------------------------------------------------------------------------
double Carac::get_image_residual(Camera* cam) const
{
  return fabs(cam->info.pix_grey - p[0] * cam->info.light_vector*n - a[0]);
}
//----------------------------------------------------------------------------*/


// Color Carac
/*Carac::Carac(Camera **cam)
{ 
  Matrix<double> A(4,4);
	Vector<double> B(4);
	
  Matrix<double> C(12,6);
	Vector<double> D(12);	
       
  for(int i=0; i<4; i++)
  {
    A(i,0) = cam[i]->info.light_vector[0];
    A(i,1) = cam[i]->info.light_vector[1];
    A(i,2) = cam[i]->info.light_vector[2];
    A(i,3) = 1.;  
    B[i]   = cam[i]->info.pix_grey/LIGHT_SOURCE;
    C(i*3+0,3) = 1.; 
    C(i*3+1,4) = 1.;
    C(i*3+2,5) = 1.;        
    D[i*3+0] = cam[i]->info.pix_color[0]/LIGHT_SOURCE_R;
    D[i*3+1] = cam[i]->info.pix_color[1]/LIGHT_SOURCE_G;
    D[i*3+2] = cam[i]->info.pix_color[2]/LIGHT_SOURCE_B;        
  }

  A = inverse(A);
        
  if(norm(A) == 0.)  // out of a picture or A not invertible
  {
    for(int i = 0; i<3; i++)
    {
      n[i] = 0.;
      a[i] = 0.;
      p[i] = 0.;
    }
  }
  else
  {
    B = A*B;
    n[0] = B[0];
    n[1] = B[1];
    n[2] = B[2];
    n /= norm(n);
    
    for(int i=0; i<4; i++)
    {
      double coeff = cam[i]->info.light_vector*n;
      C(i*3+0,0) = coeff;
      C(i*3+1,1) = coeff;
      C(i*3+2,2) = coeff;
    }
    
    C = pinverse(C);
   
    if(norm(C) == 0.)  // SVD fails
    {
      for(int i = 0; i<3; i++)
      {
        n[i] = 0.;
        a[i] = 0.;
        p[i] = 0.;
      }
    }
    else
    {
    	D = C * D;
    	p[0] = D[0];   	p[1] = D[1];   	p[2] = D[2];
    	a[0] = D[3];   	a[1] = D[4];   	a[2] = D[5];  	
    }
  }
}
//----------------------------------------------------------------------------
double Carac::get_image_residual(Camera* cam) const
{
  double coeff = cam->info.light_vector*n;
  double r = 0.;  
  for(int c=0; c<3; c++)
    r += fabs(cam->info.pix_color[c] - p[c] * coeff - a[c]);
  return r;
}
//----------------------------------------------------------------------------*/


// p, a -> n
/*Carac::Carac(Camera **cam)
{ 
  Matrix<double> A(12,12);
	Vector<double> B(12);
  
  for(int i=0; i<12; i++)
    for(int j=0; j<12; j++)  
      A(i,j) = 0.;
       
  for(int i=0; i<4; i++)
  {
    for(int c=0; c<3; c++)
    {
      A(i*3+c,c*3+0) = cam[i]->info.light_vector[0];
      A(i*3+c,c*3+1) = cam[i]->info.light_vector[1];
      A(i*3+c,c*3+2) = cam[i]->info.light_vector[2];
      A(i*3+c,9+c) = 1.;
      B[i*3+c] = cam[i]->info.pix_color[c]/(c==0?LIGHT_SOURCE_R:(c==1?LIGHT_SOURCE_G:LIGHT_SOURCE_B));  
    }
  }

  A = inverse(A);
        
  if(norm(A) == 0.)  // out of a picture or A not invertible
  {
    n = vec(0.,0.,0.);
    p = vec(0.,0.,0.);
    a = vec(0.,0.,0.);    
  }
  else
  {   
    B = A*B;
    
    for(int c=0; c<3; c++)
    {
      vec nc = vec(B[c*3+0], B[c*3+1], B[c*3+2]);
      p[c] = norm(nc);
      a[c] = B[9+c]*(c==0?LIGHT_SOURCE_R:(c==1?LIGHT_SOURCE_G:LIGHT_SOURCE_B));  
    }
    
    Matrix<double> C(12,3);

    for(int i=0; i<4; i++)
    {
      for(int c=0; c<3; c++)
      {
        vec l = cam[i]->info.light_vector*p[c]*(c==0?LIGHT_SOURCE_R:(c==1?LIGHT_SOURCE_G:LIGHT_SOURCE_B));
        C(i*3+c,0) = l[0];
        C(i*3+c,1) = l[1];
        C(i*3+c,2) = l[2];
        B[i*3+c] = cam[i]->info.pix_color[c] - a[c];  
      }
    }

    C = pinverse(C);
   
    if(norm(C) == 0.)  // SVD fails
    {
      for(int i = 0; i<3; i++)
      {
        n[i] = 0.;
        a[i] = 0.;
        p[i] = 0.;
      }
    }
    else
    {
    	Vector<double> N = C * B;
    	n[0] = N[0];
    	n[1] = N[1];
    	n[2] = N[2];    	    	
    }
  }
}
//----------------------------------------------------------------------------
double Carac::get_image_residual(Camera* cam) const
{
  double coeff = cam->info.light_vector * n;
  double r = 0.;  
  for(int c=0; c<3; c++)
    r += fabs(cam->info.pix_color[c] - (c==0?LIGHT_SOURCE_R:(c==1?LIGHT_SOURCE_G:LIGHT_SOURCE_B)) * p[c] * coeff - a[c]);
  return r;
}
//----------------------------------------------------------------------------*/





