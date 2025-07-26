#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <stdbool.h>

typedef void (*ODE)(double,gsl_vector*,gsl_vector*);
typedef bool (*Condizione)(double,gsl_vector*);

struct InfoBaseSimulazione{
  ODE dinamica;
  Condizione condizione;
  double *tCondizione;
  size_t *indiceCondizione;
  double t0,T,h;
};

gsl_matrix* EuleroAvanti(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* EuleroIndietro(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* CrankNicolson(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* Heun(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* RungeKuttaEsplicito(struct InfoBaseSimulazione* infoSimulazione,double* A_Butcher, double* B_Butcher,const unsigned stadi,gsl_vector* statoIniziale);
gsl_matrix* LMM(struct InfoBaseSimulazione* infoSimulazione,double* A_LMM,double* B_LMM,double b_1,gsl_matrix* innesco);
