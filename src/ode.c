#include <math.h>
#include <gsl/gsl_blas.h>
#include "../include/ode.h"

gsl_matrix* EuleroAvanti(void (*f_ODE)(double,gsl_vector*,gsl_vector*),const double T,const double h, const gsl_vector* y0){
  const unsigned N_h=(unsigned)(floor(T/h)+1.0);
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(y0->size,N_h);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,y0);
  
  //EA
  gsl_vector* dy_Buffer=gsl_vector_alloc(y0->size);
  double t=h;
  for(unsigned k=1;k < N_h; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_add(&(O_k.vector), &(O_k_1.vector));
    
    f_ODE(t-h,&(O_k_1.vector),dy_Buffer);
    gsl_vector_scale( dy_Buffer,h );
    gsl_vector_add( &(O_k.vector), dy_Buffer);
    t += h;
  }
  gsl_vector_free(dy_Buffer);
  return O_sim;
}

gsl_matrix* EuleroIndietro(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0){
  const unsigned N_h=(unsigned)(floor(T/h)+1.0);
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(y0->size,N_h);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,y0);
  
  //EI
  gsl_vector* dy_Buffer=gsl_vector_alloc(y0->size);
  gsl_vector* PF_tmp=gsl_vector_alloc(y0->size);
  double t=h;
  for(unsigned k=1;k < N_h; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    
    //Iterazioni di punto fisso
    double epsilonPF=1e-14,errore=GSL_POSINF;
    unsigned MAX_ITERS=50;
    unsigned k=1;
    while(errore >= epsilonPF && k <= MAX_ITERS){    
      f_ODE(t,&(O_k.vector),dy_Buffer);
      gsl_vector_scale( dy_Buffer,h );
      gsl_vector_add(dy_Buffer,&(O_k_1.vector));
      gsl_vector_memcpy(PF_tmp,&(O_k.vector));
      gsl_vector_memcpy(&(O_k.vector),dy_Buffer);
      
      //Uso by_Buffer per il calcolo dell'errore
      gsl_vector_sub(PF_tmp,dy_Buffer);
      errore=gsl_blas_dnrm2(PF_tmp);
      ++k;
    }
    t += h;
  }
  gsl_vector_free(dy_Buffer);
  gsl_vector_free(PF_tmp);
  return O_sim;
}

gsl_matrix* CrankNicolson(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0){
  const unsigned N_h=(unsigned)(floor(T/h)+1.0);
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(y0->size,N_h);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,y0);
  
  //CN
  gsl_matrix* dy_Buffer=gsl_matrix_alloc(y0->size,2);
  gsl_vector* PF_tmp=gsl_vector_alloc(y0->size);
  double t=h;
  for(unsigned k=1;k < N_h; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    
    gsl_vector_view f_k_1=gsl_matrix_column(dy_Buffer,0);
    gsl_vector_view f_k=gsl_matrix_column(dy_Buffer,1);
    f_ODE(t-h,&(O_k_1.vector),&(f_k_1.vector));
    
    //Iterazioni di punto fisso
    double epsilonPF=1e-14,errore=GSL_POSINF;
    unsigned MAX_ITERS=50;
    unsigned k=1;
    while(errore >= epsilonPF && k <= MAX_ITERS){    
      f_ODE(t,&(O_k.vector),&(f_k.vector));
      gsl_vector_add(&(f_k.vector),&(f_k_1.vector));
      
      gsl_vector_scale( &(f_k.vector),h/2 );
      gsl_vector_add( &(f_k.vector),&(O_k_1.vector));
      gsl_vector_memcpy(PF_tmp,&(O_k.vector));
      gsl_vector_memcpy(&(O_k.vector),&(f_k.vector));
      
      //Uso by_Buffer per il calcolo dell'errore
      gsl_vector_sub(PF_tmp,&(f_k.vector));
      errore=gsl_blas_dnrm2(PF_tmp);
      ++k;
    }
    t += h;
  }
  gsl_matrix_free(dy_Buffer);
  gsl_vector_free(PF_tmp);
  return O_sim;
}

gsl_matrix* Heun(void (*f_ODE)(double,gsl_vector*,gsl_vector*),const double T,const double h, const gsl_vector* y0){
  const unsigned N_h=(unsigned)(floor(T/h)+1.0);
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(y0->size,N_h);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,y0);
  
  //Heun
  gsl_matrix* dy_Buffer=gsl_matrix_alloc(y0->size,2);
  double t=h;
  for(unsigned k=1;k < N_h; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    
    gsl_vector_view f_k_1=gsl_matrix_column(dy_Buffer,0);
    gsl_vector_view f_k=gsl_matrix_column(dy_Buffer,1);
    f_ODE(t-h,&(O_k_1.vector),&(f_k_1.vector));
    
    
    //Calcolo EA
    gsl_vector_scale( &(f_k_1.vector),h );
    gsl_vector_add( &(O_k.vector), &(f_k_1.vector));
    f_ODE(t,&(O_k.vector),&(f_k.vector));
    
    //Calcolo passo successivo
    gsl_vector_scale( &(f_k.vector),h );
    gsl_vector_add( &(f_k.vector), &(f_k_1.vector));
    gsl_vector_scale( &(f_k.vector),5e-1 );
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    gsl_vector_add( &(O_k.vector), &(f_k.vector));
    t += h;
  }
  gsl_matrix_free(dy_Buffer);
  return O_sim;
}

gsl_matrix* RungeKuttaEsplicito(void (*f_ODE)(double,gsl_vector*,gsl_vector*),double T, double h, const gsl_vector* y0, gsl_matrix* A, gsl_vector* B){
  const unsigned N_h=(unsigned)(floor(T/h)+1.0);
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(y0->size,N_h);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,y0);
  
  //RungeKutta
  const unsigned stadi=B->size;
  gsl_matrix* K=gsl_matrix_calloc(y0->size,stadi);
  gsl_vector* C=gsl_vector_calloc(stadi);
  gsl_vector* f_k=gsl_vector_calloc(y0->size);
  gsl_vector* uni=gsl_vector_calloc(stadi);
  gsl_vector_set_all(uni,1.0);
  
  //Ottengo i coefficienti c_j
  gsl_blas_dgemv(CblasNoTrans,1.0,A,uni,0.0,C);
  
  double t=h;
  for(unsigned k=1;k < N_h; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    gsl_matrix_set_zero(K);
    
    //Calcolo dei K
    for(unsigned j=0; j<stadi; ++j){
      gsl_vector_view a_j=gsl_matrix_row(A,j);
      double t_j=t+h*(gsl_vector_get(C,j)-1);
      
      gsl_blas_dgemv(CblasNoTrans,1.0,K,&(a_j.vector),0.0,f_k);
      gsl_vector_scale( f_k,h );
      gsl_vector_add( &(O_k.vector), f_k);
      f_ODE(t_j,&(O_k.vector),f_k);
      gsl_matrix_set_col(K,j,f_k);
      
      gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    }
    
    //Calcolo passo successivo
    gsl_blas_dgemv(CblasNoTrans,1.0,K,B,0.0,f_k);
    gsl_vector_scale(f_k,h);
    gsl_vector_add( &(O_k.vector),f_k);
    t += h;
  }
  gsl_matrix_free(K);
  gsl_vector_free(uni);
  gsl_vector_free(C);
  gsl_vector_free(f_k);
  return O_sim;
}
