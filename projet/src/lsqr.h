#ifndef LSQRH
#define LSQRH

#include <map>
#include <utility>
#include "carac.h"

class Zvect
{
  private:
  void destroy();
  void normalize();
  
  protected:
  int _size;
  long double *data;
  int *n_shared;
  
  public:
  Zvect(int n, long double *d = NULL);
  Zvect(const Zvect& Z);
  ~Zvect();
  
  void operator = (const Zvect& Z);
  Zvect operator *(const long double& s) const; 
  friend Zvect operator *(const long double& s, const Zvect &Z);
  Zvect operator + (const Zvect& Z);
  Zvect operator - (const Zvect& Z);
    
  inline long double operator[] (int i) const {return norm*data[i];}
 
  long double norm; 
  inline int size() const {return _size;}
};



class SparseMat
{
  private:
  typedef map<int, long double> rowcol;
  vector<rowcol> rows;
  vector<rowcol> cols;  
  long double sum_norm;
  
  protected:
  int n, m;
  
  public:
  SparseMat(int n_rows, int n_cols);
  void add_coeff(int row, int col, long double v);
  void add_row_coeff(int row, list<pair<int, long double> > coeff);

  inline int num_row() {return n;}
  inline int num_col() {return m;}   
  
  class ZvectProd
  {
    private:
    const int n, m;
    vector<rowcol> &data;
    
    public:
    ZvectProd(int n_row, int n_col, vector<rowcol> &d):n(n_row),m(n_col),data(d) {}
    Zvect operator * (const Zvect& Z);
    inline int num_row() {return n;}
    inline int num_col() {return m;}    
  };
  
  ZvectProd get_prod_iter(bool transpose);
  inline long double norm();
};

typedef pair<SparseMat*, Zvect*> System;

class Surface
{
  private:
  int w, h;
  const long double fx, fy;
  Image<long double>& depth_map;
  Image<Carac>& caracs;  
  Image<long double> mu;
  
  list<pair<int, long double> > make_filter(int u, int v, long double *filter,const long double& weight);

  public:
  Surface(Image<long double>& _depth_map, long double fx, long double fy, Image<Carac>& _caracs);
  
  System allocate_system();
};

Zvect solve_system(bool verbose, System S, long double BTOL = 0.01, long double ATOL = 0.);

#endif
