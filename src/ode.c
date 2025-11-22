#include <math.h>
#include <gsl/gsl_blas.h>
#include <ode.h> 

/*! \brief Struttura dati per le informazioni sulla matrice del calcolo
 *
 *  Questa matrice viene scritta sempre all'inizio del file binario per specificare la disposizione del dati nel file
 */
struct InfoMatrice{
  size_t dimensioneStruct; /*!< Dimensione di questa struttura*/
  size_t righe; /*!< Numero di righe della matrice, ovvero la dimensione dello stato*/
  size_t colonne; /*!< Numero di colonne della matrice, ovvero numero di passi del calcolo*/
  double h; /*!< Passo di integrazione*/
  double T; /*!< Periodo di integrazione*/
  double t0; /*!< Istante iniziale del calcolo*/
};

gsl_matrix* EuleroAvanti(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale){
  const size_t NumeroCampioni=(size_t)floor(infoSimulazione->T/infoSimulazione->h)+1;
  const size_t n=statoIniziale->size;
  //Allocazione matrice e inserimento stato iniziale
  gsl_matrix* O_sim=gsl_matrix_alloc(n,NumeroCampioni);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,statoIniziale);
  
  //EA
  double f_Buffer[n];
  double t_k_1=infoSimulazione->t0;
  gsl_vector_view dy_Buffer=gsl_vector_view_array(f_Buffer,n);
  size_t k=1;
  for(;k < NumeroCampioni; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    //Se è verificata la condizione termino
    bool verificaCondizione = infoSimulazione->condizione == NULL ? 0 : infoSimulazione->condizione(t_k_1,&(O_k_1.vector));
    if(verificaCondizione) break;
    
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_add(&(O_k.vector), &(O_k_1.vector));
    
    infoSimulazione->dinamica(t_k_1,&(O_k_1.vector),&(dy_Buffer.vector));
    gsl_vector_scale( &(dy_Buffer.vector),infoSimulazione->h );
    gsl_vector_add( &(O_k.vector), &(dy_Buffer.vector));
    t_k_1 = infoSimulazione->t0+((double)k)*infoSimulazione->h;
  }
  
  //Assegno gli istanti della condizione nel caso sia verificata
  if(infoSimulazione->tCondizione) *(infoSimulazione->tCondizione)=t_k_1;
  if(infoSimulazione->indiceCondizione) *(infoSimulazione->indiceCondizione)=k-1;
  return O_sim;
}

gsl_matrix* EuleroIndietro(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale){
  const size_t NumeroCampioni=(size_t)floor(infoSimulazione->T/infoSimulazione->h)+1;
  const size_t n=statoIniziale->size;
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(n,NumeroCampioni);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,statoIniziale);
  
  //EI
  double f_Buffer[n], puntoFisso_tmp[n];
  gsl_vector_view dy_Buffer=gsl_vector_view_array(f_Buffer,n);
  gsl_vector_view PF_tmp=gsl_vector_view_array(puntoFisso_tmp,n);
  double t_k=infoSimulazione->t0+infoSimulazione->h, t_k_1=infoSimulazione->t0;
  unsigned k=1;
  for(;k < NumeroCampioni; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    //Se è verificata la condizione termino
    bool verificaCondizione = infoSimulazione->condizione == NULL ? 0 : infoSimulazione->condizione(t_k_1,&(O_k_1.vector));
    if(verificaCondizione) break;
    
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    
    //Iterazioni di punto fisso
    double epsilonPF=1e-15,errore=GSL_POSINF;
    unsigned MAX_ITERS=50;
    unsigned j=1;
    while(errore >= epsilonPF && j <= MAX_ITERS){    
      infoSimulazione->dinamica(t_k,&(O_k.vector),&(dy_Buffer.vector));
      gsl_vector_scale(&(dy_Buffer.vector),infoSimulazione->h);
      gsl_vector_add(&(dy_Buffer.vector),&(O_k_1.vector));
      gsl_vector_memcpy(&(PF_tmp.vector),&(O_k.vector));
      gsl_vector_memcpy(&(O_k.vector),&(dy_Buffer.vector));
      
      //Uso by_Buffer per il calcolo dell'errore
      gsl_vector_sub(&(PF_tmp.vector),&(dy_Buffer.vector));
      errore=gsl_blas_dnrm2(&(PF_tmp.vector));
      ++j;
    }
    t_k_1=t_k;
    t_k = infoSimulazione->t0+((double)(k+1))*infoSimulazione->h;
  }
  //Assegno gli istanti della condizione nel caso sia verificata
  if(infoSimulazione->tCondizione) *(infoSimulazione->tCondizione)=t_k_1;
  if(infoSimulazione->indiceCondizione) *(infoSimulazione->indiceCondizione)=k-1;
  return O_sim;
}

gsl_matrix* CrankNicolson(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale){
  const size_t NumeroCampioni=(size_t)floor(infoSimulazione->T/infoSimulazione->h)+1;
  const size_t n=statoIniziale->size;
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(n,NumeroCampioni);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,statoIniziale);
  
  //CN
  double f_Buffer[n*2], puntoFisso_tmp[n];
  gsl_matrix_view dy_Buffer=gsl_matrix_view_array(f_Buffer,n,2);
  gsl_vector_view PF_tmp=gsl_vector_view_array(puntoFisso_tmp,n);
  double t_k=infoSimulazione->t0+infoSimulazione->h,t_k_1=infoSimulazione->t0;
  unsigned k=1;
  for(;k < NumeroCampioni; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);    
    //Se è verificata la condizione termino
    bool verificaCondizione = infoSimulazione->condizione == NULL ? 0 : infoSimulazione->condizione(t_k_1,&(O_k_1.vector));
    if(verificaCondizione) break;
    
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    gsl_vector_view f_k_1=gsl_matrix_column(&(dy_Buffer.matrix),0);
    gsl_vector_view f_k=gsl_matrix_column(&(dy_Buffer.matrix),1);
    infoSimulazione->dinamica(t_k_1,&(O_k_1.vector),&(f_k_1.vector));
    
    //Iterazioni di punto fisso
    double epsilonPF=1e-15,errore=GSL_POSINF;
    unsigned MAX_ITERS=50;
    unsigned j=1;
    while(errore >= epsilonPF && j <= MAX_ITERS){    
      infoSimulazione->dinamica(t_k,&(O_k.vector),&(f_k.vector));
      gsl_vector_add(&(f_k.vector),&(f_k_1.vector));
      
      gsl_vector_scale( &(f_k.vector),(infoSimulazione->h)/2.0 );
      gsl_vector_add( &(f_k.vector),&(O_k_1.vector));
      gsl_vector_memcpy(&(PF_tmp.vector),&(O_k.vector));
      gsl_vector_memcpy(&(O_k.vector),&(f_k.vector));
      
      //Uso by_Buffer per il calcolo dell'errore
      gsl_vector_sub(&(PF_tmp.vector),&(f_k.vector));
      errore=gsl_blas_dnrm2(&(PF_tmp.vector));
      ++j;
    }
    t_k_1=t_k;
    t_k = infoSimulazione->t0+((double)(k+1))*infoSimulazione->h;
  }
  //Assegno gli istanti della condizione nel caso sia verificata
  if(infoSimulazione->tCondizione) *(infoSimulazione->tCondizione)=t_k_1;
  if(infoSimulazione->indiceCondizione) *(infoSimulazione->indiceCondizione)=k-1;
  return O_sim;
}

gsl_matrix* Heun(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale){
  const size_t NumeroCampioni=(size_t)floor(infoSimulazione->T/infoSimulazione->h)+1;
  const size_t n=statoIniziale->size;
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(n,NumeroCampioni);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,statoIniziale);
  
  //Heun
  double f_Buffer[n*2];
  gsl_matrix_view dy_Buffer=gsl_matrix_view_array(f_Buffer,n,2);
  double t_k=infoSimulazione->t0+infoSimulazione->h,t_k_1=infoSimulazione->t0;
  unsigned k=1;
  for(;k < NumeroCampioni; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    //Se è verificata la condizione termino
    bool verificaCondizione = infoSimulazione->condizione == NULL ? 0 : infoSimulazione->condizione(t_k_1,&(O_k_1.vector));
    if(verificaCondizione) break;
    
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    gsl_vector_view f_k_1=gsl_matrix_column(&(dy_Buffer.matrix),0);
    gsl_vector_view f_k=gsl_matrix_column(&(dy_Buffer.matrix),1);
    infoSimulazione->dinamica(t_k_1,&(O_k_1.vector),&(f_k_1.vector));
    
    //Calcolo EA
    gsl_vector_scale(&(f_k_1.vector),infoSimulazione->h);
    gsl_vector_add(&(O_k.vector),&(f_k_1.vector));
    infoSimulazione->dinamica(t_k,&(O_k.vector),&(f_k.vector));
    
    //Calcolo passo successivo
    gsl_vector_scale(&(f_k.vector),infoSimulazione->h);
    gsl_vector_add(&(f_k.vector),&(f_k_1.vector));
    gsl_vector_scale(&(f_k.vector),5e-1);
    gsl_vector_memcpy(&(O_k.vector),&(O_k_1.vector));
    gsl_vector_add(&(O_k.vector),&(f_k.vector));
    t_k_1=t_k;
    t_k = infoSimulazione->t0+((double)(k+1))*infoSimulazione->h;
  }
  //Assegno gli istanti della condizione nel caso sia verificata
  if(infoSimulazione->tCondizione) *(infoSimulazione->tCondizione)=t_k_1;
  if(infoSimulazione->indiceCondizione) *(infoSimulazione->indiceCondizione)=k-1;
  return O_sim;
}

gsl_matrix* RungeKuttaEsplicito(struct InfoBaseSimulazione* infoSimulazione,double* A_Butcher, double* B_Butcher,const unsigned stadi,gsl_vector* statoIniziale){
  const size_t NumeroCampioni=(size_t)floor(infoSimulazione->T/infoSimulazione->h)+1;
  const size_t n=statoIniziale->size;
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(n,NumeroCampioni);
  gsl_matrix_set_zero(O_sim);
  gsl_matrix_set_col(O_sim,0,statoIniziale);
  
  //RungeKutta
  double K_mat[n*stadi],C_vec[stadi],f_k_vec[n],uni_vec[stadi];
  gsl_matrix_view K=gsl_matrix_view_array(K_mat,n,stadi);
  gsl_vector_view C=gsl_vector_view_array(C_vec,stadi);
  gsl_vector_view f_k=gsl_vector_view_array(f_k_vec,n);
  gsl_vector_view uni=gsl_vector_view_array(uni_vec,stadi);
  gsl_matrix_view A=gsl_matrix_view_array(A_Butcher,stadi,stadi);
  gsl_vector_view B=gsl_vector_view_array(B_Butcher,stadi);
  gsl_vector_set_all(&(uni.vector),1.0);
  
  //Ottengo i coefficienti c_j
  gsl_blas_dgemv(CblasNoTrans,1.0,&(A.matrix),&(uni.vector),0.0,&(C.vector));
  
  double t_k=infoSimulazione->t0+infoSimulazione->h,t_k_1=infoSimulazione->t0;
  unsigned k=1;
  for(;k < NumeroCampioni; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    //Se è verificata la condizione termino
    bool verificaCondizione = infoSimulazione->condizione == NULL ? 0 : infoSimulazione->condizione(t_k_1,&(O_k_1.vector));
    if(verificaCondizione) break;
    
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    gsl_matrix_set_zero(&(K.matrix));
    
    //Calcolo dei K
    for(unsigned j=0; j<stadi; ++j){
      gsl_vector_view a_j=gsl_matrix_row(&(A.matrix),j);
      double t_j=t_k_1+(infoSimulazione->h*gsl_vector_get(&(C.vector),j));
      
      gsl_blas_dgemv(CblasNoTrans,1.0,&(K.matrix),&(a_j.vector),0.0,&(f_k.vector));
      gsl_vector_scale(&(f_k.vector),infoSimulazione->h);
      gsl_vector_add(&(O_k.vector),&(f_k.vector));
      infoSimulazione->dinamica(t_j,&(O_k.vector),&(f_k.vector));
      gsl_matrix_set_col(&(K.matrix),j,&(f_k.vector));
      
      gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    }
    
    //Calcolo passo successivo
    gsl_blas_dgemv(CblasNoTrans,1.0,&(K.matrix),&(B.vector),0.0,&(f_k.vector));
    gsl_vector_scale(&(f_k.vector),infoSimulazione->h);
    gsl_vector_add(&(O_k.vector),&(f_k.vector));
    t_k_1=t_k;
    t_k = infoSimulazione->t0+((double)(k+1))*infoSimulazione->h;
  }
  //Assegno gli istanti della condizione nel caso sia verificata
  if(infoSimulazione->tCondizione) *(infoSimulazione->tCondizione)=t_k_1;
  if(infoSimulazione->indiceCondizione) *(infoSimulazione->indiceCondizione)=k-1;
  return O_sim;
}

gsl_matrix* LMM(struct InfoBaseSimulazione* infoSimulazione,double* A_LMM,double* B_LMM,double b_1,gsl_matrix* innesco){
  const size_t NumeroCampioni=(size_t)floor(infoSimulazione->T/infoSimulazione->h)+1;
  const size_t n=innesco->size1;
  const size_t p=innesco->size2-1;
  
  //Allocazione matrice
  gsl_matrix* O_sim=gsl_matrix_alloc(n,NumeroCampioni);
  gsl_matrix_set_zero(O_sim);
  
  //LMM
  //Inizializzazione dei buffer
  double Buffer_O_mat[n*(p+1)],Buffer_F_mat[n*(p+1)],CombA_vec[n],CombB_vec[n],f_1_vec[n],puntoFisso_tmp[n];
  gsl_matrix_view Buffer_O=gsl_matrix_view_array(Buffer_O_mat,n,p+1);
  gsl_matrix_view Buffer_F=gsl_matrix_view_array(Buffer_F_mat,n,p+1);
  gsl_vector_view CombA=gsl_vector_view_array(CombA_vec,n);
  gsl_vector_view CombB=gsl_vector_view_array(CombB_vec,n);
  gsl_vector_view f_1=gsl_vector_view_array(f_1_vec,n);
  gsl_vector_view PF_tmp=gsl_vector_view_array(puntoFisso_tmp,n);
  gsl_vector_view A=gsl_vector_view_array(A_LMM,p+1);
  gsl_vector_view B=gsl_vector_view_array(B_LMM,p+1);
  
  size_t k=0;
  for(; k<p+1; ++k){
    gsl_vector_view col_k_O=gsl_matrix_column(innesco,k);
    gsl_vector_view col_p_k_F=gsl_matrix_column(&(Buffer_F.matrix),p-k);
    gsl_matrix_set_col(O_sim,k,&(col_k_O.vector));
    
    double t_k=infoSimulazione->t0+((double)k)*infoSimulazione->h;
    infoSimulazione->dinamica(t_k,&(col_k_O.vector),&(col_p_k_F.vector));
    gsl_matrix_set_col(&(Buffer_O.matrix),p-k,&(col_k_O.vector));
  }
  double t_k=infoSimulazione->t0+((double)(k))*infoSimulazione->h,t_k_1=infoSimulazione->t0+((double)(k-1))*infoSimulazione->h;
  for(;k < NumeroCampioni; ++k){
    gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
    //Se è verificata la condizione termino
    //printf("\nO_k_1\n");
    //gsl_vector_fprintf(stdout,&(O_k_1.vector),"%.10lf");
    bool verificaCondizione = infoSimulazione->condizione == NULL ? 0 : infoSimulazione->condizione(t_k_1,&(O_k_1.vector));
    if(verificaCondizione) break;
    
    gsl_vector_view O_k=gsl_matrix_column(O_sim,k);
    gsl_blas_dgemv(CblasNoTrans,1.0,&(Buffer_O.matrix),&(A.vector),0.0,&(CombA.vector));
    gsl_blas_dgemv(CblasNoTrans,1.0,&(Buffer_F.matrix),&(B.vector),0.0,&(CombB.vector));
    gsl_vector_scale(&(CombB.vector),infoSimulazione->h);
    
    //Se esplicito
    if(b_1 == 0.0){
      gsl_vector_add(&(CombA.vector),&(CombB.vector));
      gsl_vector_add(&(O_k.vector),&(CombA.vector));
    }else{ //Se implicito
      //gsl_vector_view O_k_1=gsl_matrix_column(O_sim,k-1);
      gsl_vector_memcpy(&(O_k.vector), &(O_k_1.vector));
    
      //Iterazioni di punto fisso
      double epsilonPF=1e-15,errore=GSL_POSINF;
      unsigned MAX_ITERS=50;
      unsigned j=1;
      while(errore >= epsilonPF && j <= MAX_ITERS){    
        infoSimulazione->dinamica(t_k,&(O_k.vector),&(f_1.vector)); //Termine implicito
        gsl_vector_scale(&(f_1.vector),b_1*infoSimulazione->h);
        gsl_vector_add(&(f_1.vector),&(CombB.vector));
        gsl_vector_add(&(f_1.vector),&(CombA.vector));
        gsl_vector_memcpy(&(PF_tmp.vector),&(O_k.vector));
        gsl_vector_memcpy(&(O_k.vector),&(f_1.vector));
        
        //Uso f_1 per il calcolo dell'errore
        gsl_vector_sub(&(PF_tmp.vector),&(f_1.vector));
        errore=gsl_blas_dnrm2(&(PF_tmp.vector));
        ++j;
      }
    }
    
    //swap buffers
    for(unsigned j=0; j<p; ++j){
      gsl_matrix_swap_columns(&(Buffer_O.matrix),p-j,p-j-1);
      gsl_matrix_swap_columns(&(Buffer_F.matrix),p-j,p-j-1);
    }
    infoSimulazione->dinamica(t_k,&(O_k.vector),&(f_1.vector));
    gsl_matrix_set_col(&(Buffer_O.matrix),0,&(O_k.vector));
    gsl_matrix_set_col(&(Buffer_F.matrix),0,&(f_1.vector));
    
    t_k_1=t_k;
    t_k = infoSimulazione->t0+((double)(k+1))*infoSimulazione->h;
  }
  //Assegno gli istanti della condizione nel caso sia verificata
  if(infoSimulazione->tCondizione) *(infoSimulazione->tCondizione)=t_k_1;
  if(infoSimulazione->indiceCondizione) *(infoSimulazione->indiceCondizione)=k-1;
  return O_sim;
}

int fwrite_matrix(FILE* file, gsl_matrix* matrice,double h, double T, double t0){
  struct InfoMatrice infoSimulazione={
    .dimensioneStruct=sizeof(struct InfoMatrice),
    .righe=matrice->size1,
    .colonne=matrice->size2,
    .h=h,
    .T=T,
    .t0=t0
  };
  if(fwrite(&infoSimulazione,sizeof(struct InfoMatrice),1,file) == 0) return 1;
  if(gsl_matrix_fwrite(file,matrice) != 0) return 1;
  return 0;
}
