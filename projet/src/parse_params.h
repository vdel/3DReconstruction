#ifndef PARSE_PARAMSH
#define PARSE_PARAMSH

/****** Light source parameters ******/
extern double LIGHT_SOURCE; // Light source intensity at unit distance (default = 1.)
extern double LIGHT_SOURCE_R;
extern double LIGHT_SOURCE_G;      
extern double LIGHT_SOURCE_B;

extern double LIGHT_X;     // Light source x coordinate  (default = 0.)
extern double LIGHT_Y;     // Light source y coordinate  (default = 0.)
extern double LIGHT_Z;     // Light source z coordinate  (default = 0.)

/****** Ep parameters ******/
extern double OutLierRatio;  // Outlier ratio for ransac         (default = 0.5)
extern double MinProba;      // MinProba for ransac              (default = 0.7)
extern double TAU;           // Residual threshold for ransac    (default = [6.0,8.0])
extern double ETA;           // Normalisation factor             (default = 1./TAU)
extern double C0;            // Energy value when ransac fails   (default = 1.)

/****** En parameters ******/
extern double Tj;            // Maximum depth difference allowed (default = 3.)
extern double C1;            // Energy value when failure        (default = 5.)

/****** Energy parameters ******/
extern double LAMBDA_N;      // Scale for En                     (default = 7.5)  
extern double LAMBDA_S;      // Scale for Es                     (default = [1.5,3.0])          

/****** Optimisation parameters ******/
extern double LAMBDA_1;      // Scale for attach to data         (default = [0.01,0.1])
extern double LAMBDA_2;      //Scale for laplacian               (default = [0.5,1.5])



#include <Imagine/Images.h>
#include<list>
using namespace std;
using namespace Imagine;

//------------------------------------------------------------------------------
class val_param_t
{
  public:
  char *name;
  char *value;
  
  val_param_t(char* name_deb, int name_len, char* val_deb, char** val_fin);
  val_param_t(const val_param_t &s);
  ~val_param_t();
  
  void operator = (const val_param_t &s);
};
//------------------------------------------------------------------------------
class params_t
{
  private:
  char *FilePath;
  bool parsed;
    
  virtual void init_params();

  public:
  params_t();
  virtual ~params_t();
  list<val_param_t> param_list;
  bool   is_value_for(const char *name);
  int    int_param_value(const char *name);
  float  float_param_value(const char *name);
  double double_param_value(const char *name);    
  char*  string_param_value(const char *name);      
  bool   load_params(const char *FP);
};
//------------------------------------------------------------------------------

void parse_ini(const string& dir);
#endif
