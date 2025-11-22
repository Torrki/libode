/*! \file ode.h
 *  \brief File header per la dichiarazione delle funzioni dei metodi
 */
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <stdbool.h>
#include <stdio.h>

typedef void (*ODE)(double,gsl_vector*,gsl_vector*); /*!< Tipo di dato per le funzioni ODE, cioè la dinamica. I parametri sono in ordine tempo,stato,calcolo della derivata*/
typedef bool (*Condizione)(double,gsl_vector*); /*!< Tipo di dato per le funzioni delle condizioni di uscita. I parametri sono in ordine tempo,stato*/

/*! \brief Struttura dati per impostare il calcolo della soluzione numerica
 */
struct InfoBaseSimulazione{
  ODE dinamica; /*!< Dinamica da integrare*/
  Condizione condizione; /*!< Condizione di uscita, nullptr per non specificarla*/
  double *tCondizione; /*!< Indirizzo a una variabile per ottenere l'istante di uscita, nullptr per non specificarla*/
  size_t *indiceCondizione; /*!< Indirizzo a una variabile per ottenere l'indice di uscita, nullptr per non specificarla */
  double t0 /*!< Istante iniziale del calcolo*/
  double T /*!< Periodo di integrazione*/
  double h; /*!< Passo di integrazione*/
};

/*! \fn gsl_matrix* EuleroAvanti(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale)
 *  \brief Metodo di integrazione Eulero Avanti
 *
 *  \param infoSimulazione Indirizzo alla struttura dati per impostare il calcolo
 *  \param statoIniziale Vettore per lo stato iniziale del calcolo
 *  \return La funzione alloca una matrice gsl_matrix che contiene il risultato numerico, la matrice deve essere deallocata dall'utente
 */
 
/*! \fn gsl_matrix* EuleroIndietro(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale)
 *  \brief Metodo di integrazione Eulero Indietro
 *
 *  \param infoSimulazione Indirizzo alla struttura dati per impostare il calcolo
 *  \param statoIniziale Vettore per lo stato iniziale del calcolo
 *  \return La funzione alloca una matrice gsl_matrix che contiene il risultato numerico, la matrice deve essere deallocata dall'utente
 */
 
/*! \fn gsl_matrix* CrankNicolson(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale)
 *  \brief Metodo di integrazione Crank Nicholson
 *
 *  \param infoSimulazione Indirizzo alla struttura dati per impostare il calcolo
 *  \param statoIniziale Vettore per lo stato iniziale del calcolo
 *  \return La funzione alloca una matrice gsl_matrix che contiene il risultato numerico, la matrice deve essere deallocata dall'utente
 */
 
/*! \fn gsl_matrix* Heun(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale)
 *  \brief Metodo di integrazione Heun
 *
 *  \param infoSimulazione Indirizzo alla struttura dati per impostare il calcolo
 *  \param statoIniziale Vettore per lo stato iniziale del calcolo
 *  \return La funzione alloca una matrice gsl_matrix che contiene il risultato numerico, la matrice deve essere deallocata dall'utente
 */
 
/*! \fn gsl_matrix* RungeKuttaEsplicito(struct InfoBaseSimulazione* infoSimulazione,double* A_Butcher, double* B_Butcher,const unsigned stadi,gsl_vector* statoIniziale)
 *  \brief Metodo di integrazione Runge Kutta
 *
 *  \param infoSimulazione Indirizzo alla struttura dati per impostare il calcolo
 *  \param A_Butcher Array dei coefficienti di Butcher, devono essere memorizzate in ordine prima le righe
 *  \param B_Butcher Array dei pesi di interpolazione
 *  \param stadi Numero di stadi del metodo
 *  \param statoIniziale Vettore per lo stato iniziale del calcolo
 *  \return La funzione alloca una matrice gsl_matrix che contiene il risultato numerico, la matrice deve essere deallocata dall'utente
 */
 
/*! \fn gsl_matrix* LMM(struct InfoBaseSimulazione* infoSimulazione,double* A_LMM,double* B_LMM,double b_1,gsl_matrix* innesco)
 *  \brief Metodo di integrazione LMM
 *
 *  \param infoSimulazione Indirizzo alla struttura dati per impostare il calcolo
 *  \param A_LMM Coefficienti sulla soluzione numerica
 *  \param B_LMM Coefficienti sulla dinamica
 *  \param b_1 Coefficiente per metodi impliciti
 *  \param innesco Matrice per innesco del metodo
 *  \return La funzione alloca una matrice gsl_matrix che contiene il risultato numerico, la matrice deve essere deallocata dall'utente
 */
 
/*! \fn int fwrite_matrix(FILE* file, gsl_matrix* matrice,double h, double T, double t0)
 *  \brief Funzione per scrivere in formato binario il risultato del calcolo
 *
 *  \param file File nel quale scrivere il risultato
 *  \param matrice Matrice che contiene il risultato numerico
 *  \param h Passo di integrazione
 *  \param T Periodo di integrazione
 *  \param t0 Istante iniziale del calcolo
 */
gsl_matrix* EuleroAvanti(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* EuleroIndietro(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* CrankNicolson(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* Heun(struct InfoBaseSimulazione* infoSimulazione,gsl_vector* statoIniziale);
gsl_matrix* RungeKuttaEsplicito(struct InfoBaseSimulazione* infoSimulazione,double* A_Butcher, double* B_Butcher,const unsigned stadi,gsl_vector* statoIniziale);
gsl_matrix* LMM(struct InfoBaseSimulazione* infoSimulazione,double* A_LMM,double* B_LMM,double b_1,gsl_matrix* innesco);
int fwrite_matrix(FILE* file, gsl_matrix* matrice,double h, double T, double t0);
