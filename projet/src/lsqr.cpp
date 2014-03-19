#include "lsqr.h"
#include "parse_params.h"

/*----------------------------------------------------------------------------*/
Zvect::Zvect(int n, long double *d)
{
  n_shared = new int;
  *n_shared = 1;  
  _size = n;
  if(d)
  {
    data = d;
    norm = 1.;
  }
  else
  {
    data = new long double[_size];
    for(int i=0; i<_size; i++)
      data[i] = 0.;
    norm = 0.; 
  }
  normalize();
}
/*----------------------------------------------------------------------------*/
Zvect::Zvect(const Zvect& Z)
{
  n_shared = Z.n_shared;
  (*n_shared)++;
  _size = Z._size;
  data = Z.data;
  norm = Z.norm;    
}
/*----------------------------------------------------------------------------*/
Zvect::~Zvect()
{
  destroy();
}
/*----------------------------------------------------------------------------*/
void Zvect::operator = (const Zvect& Z)
{
  destroy();
  n_shared = Z.n_shared;
  (*n_shared)++;
  _size = Z._size;
  data = Z.data;
  norm = Z.norm;    
}
/*----------------------------------------------------------------------------*/
void Zvect::destroy()
{
  (*n_shared)--;
  if(*n_shared == 0)
  {
    delete[] data;
    delete n_shared;
  }
}
/*----------------------------------------------------------------------------*/
void Zvect::normalize()
{
  if(norm == 0.) return;
  long double data_norm = 0;
  for(int i=0; i<_size; i++)
    data_norm += data[i]*data[i];
  data_norm = sqrtl(data_norm);
  norm *= data_norm;
  
  for(int i=0; i<_size; i++)
    data[i] /= data_norm;
}
/*----------------------------------------------------------------------------*/
Zvect Zvect::operator *(const long double& s) const
{
  Zvect res(*this);
  res.norm *= s;
  return res;
}
/*----------------------------------------------------------------------------*/
Zvect operator *(const long double& s, const Zvect &Z)
{
  return Z*s;
}
/*----------------------------------------------------------------------------*/
Zvect Zvect::operator + (const Zvect& Z)
{
  assert(_size == Z._size);
  long double *d = new long double[_size];
  for(int i=0; i<_size; i++)
    d[i] = (*this)[i] + Z[i];
  return Zvect(_size, d);
}
/*----------------------------------------------------------------------------*/
Zvect Zvect::operator -(const Zvect& Z)  
{
  assert(_size == Z._size);
  long double *d = new long double[_size];
  for(int i=0; i<_size; i++)
    d[i] = (*this)[i] - Z[i];  
  return Zvect(_size, d);
}
/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
SparseMat::SparseMat(int n_rows, int n_cols)
{
  n = n_rows;
  m = n_cols;
  sum_norm = 0.;
  rows.resize(n);
  cols.resize(m);  
}
/*----------------------------------------------------------------------------*/
void SparseMat::add_coeff(int row, int col, long double v)
{
  assert(row>=0 && row<n);
  assert(col>=0 && col<m);
  
  rowcol::iterator it = rows[row].find(col);
  if(it != rows[row].end())
  {
    sum_norm -= it->second*it->second;
    v += it->second;
    sum_norm += v*v;
  }

  rows[row][col] = v;
  cols[col][row] = v;
}
/*----------------------------------------------------------------------------*/
void SparseMat::add_row_coeff(int row, list<pair<int, long double> > coeff)
{
  list<pair<int, long double> >::iterator it;
  for(it = coeff.begin(); it!=coeff.end(); it++)
    add_coeff(row, it->first, it->second);
}
/*----------------------------------------------------------------------------*/
long double SparseMat::norm()
{
  return sqrtl(sum_norm);
}
/*----------------------------------------------------------------------------*/
SparseMat::ZvectProd SparseMat::get_prod_iter(bool transpose)
{
  if(!transpose)
    return ZvectProd(n,m,rows);
  else
    return ZvectProd(m,n,cols);  
}
/*----------------------------------------------------------------------------*/
Zvect SparseMat::ZvectProd::operator * (const Zvect& Z)
{
  assert(m == Z.size());
  long double *d = new long double[n];
  
  rowcol::iterator iter;
  for(int i=0; i<n; i++)
  {
    d[i] = 0.;
    for(iter=data[i].begin(); iter!=data[i].end(); iter++)
      d[i] += Z[iter->first] * iter->second;
  }
  return Zvect(n, d);
}
/*----------------------------------------------------------------------------*/


  


/*----------------------------------------------------------------------------*/
Surface::Surface(Image<long double>& _depth_map, long double _fx, long double _fy, Image<Carac>& _caracs):
   fx(_fx),
   fy(_fy),
   depth_map(_depth_map),
   caracs(_caracs) 
{
  w = depth_map.width();
  h = depth_map.height();
  mu = Image<long double>(w,h);
    
  const long double w2 = (long double)(w)/2.;
  const long double h2 = (long double)(h)/2.;
  for(int u=0; u<w; u++)
    for(int v=0; v<h; v++)
    {
      long double x = ((long double)(u)-w2)/fx;
      long double y = ((long double)(v)-h2)/fy;
      mu(u,v) = sqrtl(x*x + y*y + 1.);
    }
}
/*----------------------------------------------------------------------------*/
System Surface::allocate_system()
{
  const int n = 4*w*h;
  const int m = w*h;
  
  const long double w2 = (long double)(w)/2.;
  const long double h2 = (long double)(h)/2.;  
  const long double sqrt_L1 = sqrtl(LAMBDA_1);
  const long double sqrt_1ML1 = sqrtl(1.-LAMBDA_1);
  const long double sqrt_L2 = sqrtl(LAMBDA_2);
    
  // builds B
  long double *b = new long double[n];
  for(int u=0; u<w; u++)
    for(int v=0; v<h; v++)
    {
      int i = u * h + v;
      b[i] = sqrt_L1 * mu(u,v) * depth_map(u,v);
    }
  for(int i=m; i<n; i++)
    b[i] = 0.;
  Zvect *B = new Zvect(n, b);
    
  // builds A
  SparseMat *A = new SparseMat(n,m);
  // part 1
  for(int u=0; u<w; u++)
    for(int v=0; v<h; v++)
    {
      int i = u * h + v;
      A->add_coeff(i, i, sqrt_L1*mu(u,v));
    }

  // part 2
  long double dx_filter[] = {-1./12., 0. , 1./12.,
                             -4./12., 0. , 4./12.,
                             -1./12., 0. , 1./12.};
  for(int u=0; u<w; u++)
    for(int v=0; v<h; v++)
    {
      int j = u * h + v;
      int i = m + j;
      long double x = ((long double)(u)-w2)/fx;
      long double y = ((long double)(v)-h2)/fy;
      const vec& n = caracs(u,v).n;
      long double dx_weight = sqrt_1ML1 * (n[2] - (n[0]*x + n[1]*y));
      A->add_row_coeff(i, make_filter(u,v, dx_filter, dx_weight));
      A->add_coeff(i, j, -sqrt_1ML1 * n[0] / fx);
    }      

  // part 3
  long double dy_filter[] = { 1./12.,  4./12. ,  1./12.,
                                0.,      0. ,     0.,
                             -1./12., -4./12. , -1./12.,};
  for(int u=0; u<w; u++)
    for(int v=0; v<h; v++)
    {
      int j = u * h + v;
      int i = 2*m + j;
      long double x = ((long double)(u)-w2)/fx;
      long double y = ((long double)(v)-h2)/fy;
      const vec& n = caracs(u,v).n;
      long double dy_weight = sqrt_1ML1 * (n[2] - (n[0]*x + n[1]*y));
      A->add_row_coeff(i, make_filter(u,v, dy_filter, dy_weight));
      A->add_coeff(i, j, -sqrt_1ML1 * n[1] / fy);
    }  

  // part 4
  long double laplacian_filter[] = { -1.,  -1. ,  -1.,
                                     -1.,   8. ,  -1.,
                                     -1.,  -1. ,  -1.};
  for(int u=0; u<w; u++)
    for(int v=0; v<h; v++)
    {
      int i = 3*m + u * h + v;
      A->add_row_coeff(i, make_filter(u,v, laplacian_filter, sqrt_L2));
    }          

  return pair<SparseMat*, Zvect*>(A,B);
}
/*----------------------------------------------------------------------------*/
list<pair<int, long double> > Surface::make_filter(int u, int v, long double *filter,const long double& weight)
{
  list<pair<int, long double> > res;
  
  int k=0;
  for(int j=v-1; j<=v+1; j++)
    for(int i=u-1; i<=u+1; i++, k++)
    {
      int x = i<0?0:(i>=w?w-1:i);
      int y = j<0?0:(j>=h?h-1:j);
      res.push_back(pair<int, long double>(x*h+y, filter[k]*weight));
    }
  return res;
}
/*----------------------------------------------------------------------------*/






/*----------------------------------------------------------------------------*/
Zvect solve_system(bool verbose, System S, long double BTOL, long double ATOL)
{
  SparseMat::ZvectProd A  = S.first->get_prod_iter(false);
  SparseMat::ZvectProd AT = S.first->get_prod_iter(true);
  const Zvect &b = *S.second;
  long double ATOLbyNORMA = ATOL * S.first->norm();
  long double BTOLbyNORMB = BTOL * b.norm;
  
  // 1 - Initialize
  Zvect u(b);
  long double beta = u.norm;
  u.norm = 1.;
  
  Zvect v(AT*u);
  long double alpha = v.norm;
  v.norm = 1.;
  
  Zvect w(v);
  Zvect x(A.num_col());
  
  long double phi_dash = beta;
  long double rho_dash = alpha;
  
  long double rho, c, s, theta, phi;
  
  int iter = 0;
  cout << endl;
  do
  {
    // 2 - Continue the bidiagonalization
    u = A*v - alpha*u;   beta = u.norm;   u.norm = 1.;
    v = AT*u - beta*v;   alpha = v.norm;  v.norm = 1.;
    
    // 3 - Construct and apply next orthogonal transformation
    rho = sqrtl(rho_dash*rho_dash + beta*beta);
    c = rho_dash / rho;
    s = beta / rho;
    theta = s * alpha;
    rho_dash = - c * alpha;
    phi = c * phi_dash;
    phi_dash = s * phi_dash;
    
    // 4 - Update x, w
    x = x + (phi/rho) * w;
    w = v - (theta/rho) * w;
    
    if(false)
    {
      cout << "||u|| = " << u.norm << endl;
      cout << "||v|| = " << v.norm << endl;    
      cout << "rho = " << rho << endl;
      cout << "c = " << c << endl;
      cout << "s = " << s << endl; 
      cout << "theta = " << theta << endl;
      cout << "rho dash = " << rho_dash << endl;
      cout << "phi = " << phi << endl;
      cout << "phi_dash = " << phi_dash << endl;        
      cout << "||x|| = " << x.norm << endl;  
      cout << "||w|| = " << w.norm << endl << endl;      
    }
    
    if(verbose)
    {      
      cout << "||r|| = " << phi_dash << endl;
      cout << "BTOL * ||B|| + ATOL * ||A|| * ||x|| = " << BTOLbyNORMB + ATOLbyNORMA * x.norm << endl;
    }
    
    iter++;
    cout << "LSQR iteration n°" << iter << ":   res = " << fabsl(s-1.) << endl;    
  }
  while(phi_dash > (BTOLbyNORMB + ATOLbyNORMA * x.norm) && fabsl(s-1.) > 0.0001);
  
  return x;
}
/*----------------------------------------------------------------------------*/
