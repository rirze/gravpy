#include <math.h>

double alpha0phir(double r, double s, double a){
  if (r/s < 1.0e-4){
    return r*pow(s,a-2.0);}
  else if (a==0.0){
    return log(1.0+r*r/(s*s))/r;}
  else { 
    return 2.0/(a*r)*(pow(s*s+r*r,a/2.0) - pow(s,a)); }
}

double alpha_integrand(int n, double args[n]){
  double u = args[0];
  double r = args[1];
  double e = args[2];
  double sint = args[3];
  double cost = args[4];
  double s = args[5];
  double a = args[6];
  int ind = (int)args[7];

  double q  = 1.0-e;
  double t0 = 1.0-(1.0-q*q)*u;
  double t1 = 1.0/sqrt(t0);
  double t3 = t1/t0;
  double t5 = t3/t0;
  double t6 = (cost*cost+sint*sint/(1.0 - (1.0-q*q)*u))*r*r;
  double t7 = t6*u;

  if (ind == 0){
    double mphiu = u!=0.0 ? sqrt(t7)*alpha0phir(sqrt(t7),s,a)/u : t6*pow(s,a-2.0) ;
    return t1*mphiu;
  }
  
  double k  = pow(s*s+t7,a/2.0-1.0);
  switch(ind){
  case 1:
    return t1*k;
  case 2:
    return t3*k;
  }
  
  double kp = k*(a/2.0-1.0)/(s*s+t7);

  switch(ind){
  case 3:
    return t1*kp*u;
  case 4:
    return t5*kp*u;
  case 5:
    return t3*kp*u;
  }

  return ind;
}
