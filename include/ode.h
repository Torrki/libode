#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

gsl_matrix* EuleroAvanti(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0);
gsl_matrix* EuleroIndietro(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0);
gsl_matrix* CrankNicolson(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0);
gsl_matrix* Heun(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0);
gsl_matrix* RungeKuttaEsplicito(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0, gsl_matrix* A, gsl_vector* B);
