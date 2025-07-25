#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <stdbool.h>

gsl_matrix* EuleroAvanti(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, double h, gsl_vector* y0,bool (*condizione)(double,gsl_vector*), double *tCondizione,unsigned* indiceCondizione);
gsl_matrix* EuleroIndietro(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, double h, gsl_vector* y0,bool (*condizione)(double,gsl_vector*), double *tCondizione,unsigned* indiceCondizione);
gsl_matrix* CrankNicolson(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, double h, gsl_vector* y0,bool (*condizione)(double,gsl_vector*), double *tCondizione,unsigned* indiceCondizione);
gsl_matrix* Heun(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, double h, gsl_vector* y0,bool (*condizione)(double,gsl_vector*), double *tCondizione,unsigned* indiceCondizione);
gsl_matrix* RungeKuttaEsplicito(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, double h, gsl_vector* y0, gsl_matrix* A, gsl_vector* B,bool (*condizione)(double,gsl_vector*), double *tCondizione,unsigned* indiceCondizione);
gsl_matrix* LMM(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double t0,double T, double h, gsl_matrix* innesco, gsl_vector* A, gsl_vector* B, double b_1,bool (*condizione)(double,gsl_vector*), double *tCondizione,unsigned* indiceCondizione);
