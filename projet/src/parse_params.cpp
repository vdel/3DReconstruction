#include "parse_params.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/****** Light source parameters ******/
double LIGHT_SOURCE; // Light source intensity at unit distance (default = 1.)
double LIGHT_SOURCE_R;
double LIGHT_SOURCE_G;      
double LIGHT_SOURCE_B;

double LIGHT_X;     // Light source x coordinate  (default = 0.)
double LIGHT_Y;     // Light source y coordinate  (default = 0.)
double LIGHT_Z;     // Light source z coordinate  (default = 0.)

/****** Ep parameters ******/
double OutLierRatio;  // Outlier ratio for ransac         (default = 0.5)
double MinProba;      // MinProba for ransac              (default = 0.7)
double TAU;           // Residual threshold for ransac    (default = [6.0,8.0])
double ETA;           // Normalisation factor             (default = 1./TAU)
double C0;            // Energy value when ransac fails   (default = 1.)

/****** En parameters ******/
double Tj;            // Maximum depth difference allowed (default = 3.)
double C1;            // Energy value when failure        (default = 5.)

/****** Energy parameters ******/
double LAMBDA_N;      // Scale for En                     (default = 7.5)  
double LAMBDA_S;      // Scale for Es                     (default = [1.5,3.0])          

/****** Optimisation parameters ******/
double LAMBDA_1;      // Scale for attach to data         (default = [0.01,0.1])
double LAMBDA_2;      //Scale for laplacian               (default = [0.5,1.5])

//------------------------------------------------------------------------------
val_param_t::val_param_t(char* name_deb, int name_len, char* val_deb, char** val_fin)
{
  while(**val_fin==' ' || **val_fin==';') (*val_fin)--;
  (*val_fin)++;
  if(name_len>0 && *val_fin-val_deb>0)
  {
    name = new char[name_len+1];
    memcpy(name,name_deb,name_len);
    name[name_len]='\0';
    value = new char[*val_fin-val_deb+1];  
    memcpy(value,val_deb,*val_fin-val_deb);
    value[*val_fin-val_deb]='\0';      
  }
  else
  {
    name=NULL;
    value=NULL;
  }
}
//------------------------------------------------------------------------------ 
val_param_t::val_param_t(const val_param_t &s)
{
  if(s.name)
  {
    name = new char[strlen(s.name)+1];
    memcpy(name,s.name,strlen(s.name)+1);
  }
  if(s.value)
  {
    value = new char[strlen(s.value)+1];
    memcpy(value,s.value,strlen(s.value)+1);
  } 
}
//------------------------------------------------------------------------------  
val_param_t::~val_param_t()
{
  if(name) delete name;
  if(value) delete value;
}
//------------------------------------------------------------------------------
void val_param_t::operator = (const val_param_t &s)
{
  if(name) {delete name; name=NULL;}
  if(value) {delete value; value=NULL;}
  if(s.name)
  {
    name = new char[strlen(s.name)+1];
    memcpy(name,s.name,strlen(s.name)+1);
  }
  if(s.value)
  {
    value = new char[strlen(s.value)+1];
    memcpy(value,s.value,strlen(s.value)+1);
  }  
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
params_t::params_t()
{
  parsed = false;
  FilePath=NULL;
  param_list.clear();
}
//------------------------------------------------------------------------------
params_t::~params_t()
{
  if(FilePath) delete FilePath;
  param_list.clear();
}
//------------------------------------------------------------------------------
void params_t::init_params()
{
}
//------------------------------------------------------------------------------
bool params_t::is_value_for(const char *name)
{
  if(!parsed) return false;
  list<val_param_t> :: iterator iter;
  for(iter = param_list.begin(); iter != param_list.end(); iter++)
  {
    if(strcmp(name,iter->name)==0)
      break;
  }
  return iter != param_list.end();
}
//------------------------------------------------------------------------------
int params_t::int_param_value(const char *name)
{
  list<val_param_t> :: iterator iter;
  for(iter = param_list.begin(); iter != param_list.end(); iter++)
    if(strcmp(name,iter->name)==0)
    {
      char *deb=iter->value;
      char *fin=&iter->value[strlen(iter->value)];
      return strtol(deb,&fin,10);
    }
  if(iter == param_list.end())
  {
    fprintf(stderr, "Paramètre manquant dans %s: %s\n",FilePath,name);
    exit(EXIT_FAILURE);
  }
  return 0;  
}
//------------------------------------------------------------------------------
float params_t::float_param_value(const char *name)
{
  list<val_param_t> :: iterator iter;
  for(iter = param_list.begin(); iter != param_list.end(); iter++)
    if(strcmp(name,iter->name)==0)
    {
      char *deb=iter->value;
      char *fin=&iter->value[strlen(iter->value)];
      return strtof(deb,&fin);
    }
  if(iter == param_list.end())
  {
    fprintf(stderr, "Paramètre manquant dans %s: %s\n",FilePath,name);
    exit(EXIT_FAILURE);
  }
  return 0.;
}
//------------------------------------------------------------------------------
double params_t::double_param_value(const char *name)
{
  list<val_param_t> :: iterator iter;
  for(iter = param_list.begin(); iter != param_list.end(); iter++)
    if(strcmp(name,iter->name)==0)
    {
      char *deb=iter->value;
      char *fin=&iter->value[strlen(iter->value)];
      return strtod(deb,&fin);
    }
  if(iter == param_list.end())
  {
    fprintf(stderr, "Paramètre manquant dans %s: %s\n",FilePath,name);
    exit(EXIT_FAILURE);
  }
  return 0.;
}
//------------------------------------------------------------------------------
char* params_t::string_param_value(const char *name)
{
  char* ptr;
  list<val_param_t> :: iterator iter;
  for(iter = param_list.begin(); iter != param_list.end(); iter++)
    if(strcmp(name,iter->name)==0)
    {
      char *deb=iter->value;
      char *fin=&iter->value[strlen(iter->value)];
      
      while(*fin==' ') fin--;
      fin++;
      ptr = new char[fin-deb+1];
      memcpy(ptr,deb,fin-deb);
      ptr[fin-deb]='\0';
      return ptr;
    }
  if(iter == param_list.end())
  {
    fprintf(stderr, "Paramètre manquant dans %s: %s\n",FilePath,name);
    exit(EXIT_FAILURE);
  }    
  return NULL;  
}
//------------------------------------------------------------------------------
bool params_t::load_params(const char *FP)
{
  FilePath=new char[strlen(FP)+1];
  memcpy(FilePath,FP,strlen(FP)+1);
  FILE *File=fopen(FilePath,"rt+");
  if(File!=NULL)
  {
    list<val_param_t> :: iterator iter;
    int Len,Line;
    char Buff[1001],*PtrDeb,*Ptr;
    Ptr=Buff;
    Len=0;
    Line=1;
    while(!feof(File))
    {
      PtrDeb=Buff;
      Len+=fread(Ptr,1,1000-Len,File);
      if(feof(File)) Buff[Len]='\n';
      while((Ptr=(char*)memchr(PtrDeb,'\n',Len-(PtrDeb-Buff))))
      {
        while(*PtrDeb==' ') PtrDeb++;
        if(*PtrDeb!=';')
        {
          while(*PtrDeb==' ' || *PtrDeb=='\r' || *PtrDeb=='\t') PtrDeb++;
          if(PtrDeb!=Ptr)
          {
            char *PtrEndVar;
            char *PtrEqual=(char*)memchr(PtrDeb,'=',Ptr-PtrDeb);
            char *PtrFin;
            if(!PtrEqual)
            {
              fprintf(stderr, "Erreur de syntaxe dans %s, ligne %d\n",FilePath,Line);
              exit(EXIT_FAILURE);
            }
            PtrEndVar = PtrEqual-1;
            PtrEqual++;
            PtrFin=PtrEqual;
            while(*PtrEndVar==' ') PtrEndVar--;            
            while(*PtrFin!=';' && PtrFin<Ptr) PtrFin++;
            if(*PtrFin=='\r') PtrFin--;
            
            for(iter = param_list.begin(); iter != param_list.end(); iter++)
              if(((unsigned int)(PtrEndVar-PtrDeb+1))>=strlen(iter->name) && 
                memcmp(PtrDeb,iter->name,((unsigned int)(PtrEndVar-PtrDeb+1)))==0)
                break;              
            // On n'a pas trouvé dans la liste des paramètres inattendus                  
            // On le rajoute
            if(iter == param_list.end())
            {
              val_param_t p(PtrDeb,PtrEndVar-PtrDeb+1,PtrEqual,&PtrFin);
              if(PtrFin==PtrEqual)
              {
                fprintf(stderr, "Valeur de paramètre incorrecte dans %s, ligne %d\n",FilePath,Line);
                exit(EXIT_FAILURE);
              }
              param_list.push_front(p);
            }
            else
            {
              fprintf(stderr, "Paramètre déjà initialisé dans %s, ligne %d\n",FilePath,Line);
              exit(EXIT_FAILURE);              
            }
          }
        }
        Line++;
        PtrDeb=Ptr+1;
      }
      Len-=(PtrDeb-Buff);
      memcpy(Buff,PtrDeb,Len);
      Ptr=&Buff[Len];
    }
    fclose(File);
    init_params();
    parsed = true;
    return true;
  }
  else
    return false;
}
//------------------------------------------------------------------------------



#define N_PARAMS 17
//------------------------------------------------------------------------------
void parse_ini(const string& dir)
{
  string file = dir;
  file += "params.ini";
  params_t P;
  P.load_params(file.c_str());  

  double *params_ptr[N_PARAMS] = {
      &LIGHT_SOURCE_R,
      &LIGHT_SOURCE_G,      
      &LIGHT_SOURCE_B,      
      &LIGHT_X,
      &LIGHT_Y,
      &LIGHT_Z,
      &OutLierRatio,
      &MinProba,
      &TAU,
      &ETA,
      &C0,
      &Tj,
      &C1,
      &LAMBDA_N,
      &LAMBDA_S,
      &LAMBDA_1,
      &LAMBDA_2};  
  const char *params_str[N_PARAMS] = {
      "LIGHT_SOURCE_R",
      "LIGHT_SOURCE_G",
      "LIGHT_SOURCE_B",            
      "LIGHT_X",
      "LIGHT_Y",
      "LIGHT_Z",
      "OutLierRatio",
      "MinProba",
      "TAU",
      "ETA",
      "C0",
      "Tj",
      "C1",
      "LAMBDA_N",
      "LAMBDA_S",
      "LAMBDA_1",
      "LAMBDA_2"};
  double params_value[N_PARAMS] = {
      1.,  //"LIGHT_SOURCE_R",
      1.,  //"LIGHT_SOURCE_G",
      1.,  //"LIGHT_SOURCE_B",            
      0.,  //"LIGHT_X",
      0.,  //"LIGHT_Y",
      0.,  //"LIGHT_Z",
      0.5, //"OutLierRatio",
      0.7, //"MinProba",      
      7.,  //"TAU",
      0.1429, //"ETA",
      1., //"C0",
      3., //"Tj",
      5., //"C1",
      7.5, //"LAMBDA_N",
      6.,  //"LAMBDA_S",
      0.1, //"LAMBDA_1",
      0.5};//"LAMBDA_2"};      
  
  for(int i=0; i<N_PARAMS; i++)
    if(P.is_value_for(params_str[i]))
      *params_ptr[i] = P.double_param_value(params_str[i]);    
    else
    {
      if(strcmp("ETA",params_str[i]) == 0)
        *params_ptr[i] = 1. / *params_ptr[i-1];
      else
        *params_ptr[i] = params_value[i];    
      printf("No value for %s, assumed %s = %f\n", params_str[i], params_str[i], *params_ptr[i]);
    } 
  
  LIGHT_SOURCE = LIGHT_SOURCE_R*0.299 +
                 LIGHT_SOURCE_G*0.587 +
                 LIGHT_SOURCE_B*0.114;
}
