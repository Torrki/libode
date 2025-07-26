#include "../../include/ode.h"

void SistemaLineare(double t, gsl_vector* y, gsl_vector* dy);

int main(int argc, char* argv[]){
  struct InfoBaseSimulazione infoSimulazione={.dinamica=SistemaLineare,.condizione=NULL,.tCondizione=NULL,.indiceCondizione=NULL,.t0=0.0,.T=1.0,.h=1e-4};
  double y0[]={2.0};
  gsl_vector_view statoIniziale=gsl_vector_view_array(y0,1);
  gsl_matrix* risultatoNumerico=Heun(&infoSimulazione,&(statoIniziale.vector));
  
  //gsl_matrix_fprintf(stdout,risultatoNumerico,"%.12lf");
  FILE* fileSimulazione=fopen("dati","wr");
  gsl_matrix_fwrite(fileSimulazione,risultatoNumerico);
  fclose(fileSimulazione);
  gsl_matrix_free(risultatoNumerico);
  return 0;
}

void SistemaLineare(double t, gsl_vector* y, gsl_vector* dy){
  gsl_vector_memcpy(dy,y);
  gsl_vector_scale(dy,-5.0);
}
