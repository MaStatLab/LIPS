#include <vector> 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_randist.h>

// #include "stdafx.h"
// #include "modelsel_pf_funcs.hh"
#include <string.h>
#include "calc.cpp"


#define MAXVAR 200
#define NUMNODEVAR 2

union INDEX_TYPE_t {
  unsigned short var[MAXVAR]; 
  unsigned long long index;
};

typedef INDEX_TYPE_t INDEX_TYPE;

typedef pair<int, double> ADD_PROB_TYPE;

bool comparator( const ADD_PROB_TYPE &l, const ADD_PROB_TYPE &r) {
    return l.second > r.second;
}

class cv_state
{
public:
  INDEX_TYPE index;
  int size;
  double logBF;
};

using namespace std;

extern double rho0;
extern double *XtX;
extern double *XtY;
extern double *XtXwork;
extern double *XtYwork;
extern double *coefficients;

extern int p;
extern int method;
extern int rho_method;
extern double alpha;
extern int nobs;
extern double yty;
extern double SSY;
extern int lIterates;
extern double **models;
extern int kstep;
extern double p_mix;
extern double rho_prop;
extern int n_top_vars;
extern int particle_id;


double get_logBF(INDEX_TYPE index,int size);
inline bool is_in_model(const INDEX_TYPE& I,int i);
inline bool is_not_in_model(const INDEX_TYPE& I, int i);
void get_model_data(INDEX_TYPE index, int size);
double compute_rho(int level, int rho_method = 0);
smc::particle<cv_state> fInitialise(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<cv_state> & part, 
	     smc::rng *pRng);
void fMove(long lTime, smc::particle<cv_state> & pFrom, 
	   smc::rng *pRng);
void fMove2(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);
void fMove3(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);
void fMove4(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);
void fMove5(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);
void fMove6(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);
void fMove7(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);

double get_add_prob(cv_state *model_ptr, int i);
inline int get_i_sub(INDEX_TYPE I, int i);
void build_subGBT(const cv_state & root_model);
INDEX_TYPE make_global_index(INDEX_TYPE I_root, INDEX_TYPE I_sub, int level);
inline INDEX_TYPE make_child_index(INDEX_TYPE& I, unsigned short add_dim);
double *get_child(INDEX_TYPE& I, int i,int level);
inline INDEX_TYPE get_next_node(INDEX_TYPE& I, int n, int k);
double *get_node(INDEX_TYPE& I, int level);
unsigned int get_node_index(INDEX_TYPE& I,int level);


unsigned long long Choose(int n, int k) {
  unsigned long long c = 1;
  unsigned long long d = 1;
  for (int i = 0; i <k; i++) {
    c *= (n-i); d *= (i+1);
  }
  return c  / d;
}

inline INDEX_TYPE init_index(int n,int level) { // What is the use of n?
    INDEX_TYPE init;
    
    for (int i=0; i < MAXVAR; i++) {
      if (i < level) init.var[i] = i+1;
      else init.var[i] = 0;
    }    
    return init;
}

inline bool is_in_model(const INDEX_TYPE& I,int i) {
    int ind=0;
    int data = i+1;
    while (ind < min(MAXVAR,p)) {
        if (I.var[ind] == data) return true;
        ind++;
    }
    return false;
}
inline bool is_not_in_model(const INDEX_TYPE& I, int i) {
    return !is_in_model(I,i);
}

// void get_model_data(const cv_state & model) {
  
//   int i,j,l;
//   INDEX_TYPE index = model.index;
//   int size = model.size;

//   for (j = 0; j < size; j++) {
//     i = index.var[j]-1;     
//     XtYwork[j] = XtY[i];
//     for (l=0; l < size; l++) {
//       XtXwork[j*size + l] = XtX[i*p + index.var[l]-1];
//     }
//   }
// }

void get_model_data(INDEX_TYPE index, int size) {
  
  int i,j,l;
  
  for (j = 0; j < size; j++) {
    i = index.var[j]-1;     
    XtYwork[j] = XtY[i];
    for (l=0; l < size; l++) {
      XtXwork[j*size + l] = XtX[i*p + index.var[l]-1];
    }
  }
}

double compute_rho(int level, int rho_method) {
  
  if (level < lIterates && level < p) {
    if (level < MAXVAR) {
      if (rho_method ==0) {
	return rho0;
      }
      
      if (rho_method ==1) {
	return 1-pow(1-rho0,level);
      }

      if (rho_method ==2) { // Beta-Binomial(1,1)
	int temp = min(min(MAXVAR,p),lIterates);
	return 1.0/(temp-level+1);
      }
    }

    else { // MAXVAR reached
     
      return 1.0;
    }
  }
  
  else {
    return 1.0;
  }
}



double get_logBF(INDEX_TYPE index, int size)
{

  //  int size = model.size;
  // double *XtXwork;
  // double *XtYwork;
  double mse_m;
  double temp;
  
  XtXwork = new double[size*size];
  XtYwork = new double[size];
  coefficients = new double[size];

  get_model_data(index,size); 
 
  mse_m = yty;
  cholreg(XtXwork, XtYwork,coefficients, &mse_m, size, nobs);
  temp = 1.0 - (mse_m * (double) ( nobs - size))/SSY;

  delete [] coefficients;
  delete [] XtXwork;
  delete [] XtYwork;
  
  //cout << mse_m << "," << nobs << "," << size << "," << SSY << ": ";
  if (method == 0) { // g-prior
    //cout << "logBF_gprior(" << temp <<"," << nobs << "," << size+1 <<"," << alpha <<"): ";
    return logBF_gprior(temp, (int)nobs, size + 1, alpha);
  } else { //hyper g-prior
#ifdef RLIB
    return logBF_hyperGprior(temp, nobs, size+1, alpha);
#else
    cout << "method =1 is not implemented" << endl;
#endif
  }
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<cv_state> fInitialise(smc::rng *pRng)
{
  cv_state model;
  
  model.index = init_index(p,0); // each particle starts with the null model
  model.size = 0;
  

  return smc::particle<cv_state>(model,0);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng) // Bootstrap filter
{
  cv_state * cv_to = pFrom.GetValuePointer();
  
  if (cv_to->size < lTime - 1 || cv_to->size >= lIterates) {  
  }
  
  //  else if (pRng->UniformS() < compute_rho(cv_to->size)) {
  else if (unifRand() < compute_rho(cv_to->size,rho_method)) {
  }
    
  else {
    // double u = pRng->UniformS();
    double u = unifRand();
    double prev_logBF;
    int i;
    double incr = 1.0/(p-cv_to->size);
    double cum_prob_curr = 0;
    
 
    for(i = 0; i < p && cum_prob_curr < u; i++  ) {
 
      if (is_not_in_model(cv_to->index,i)) {
	cum_prob_curr += incr;

      }
    }

    cv_to->index = make_child_index(cv_to->index,i-1);
    cv_to->size += 1;
    prev_logBF = cv_to->logBF;

    cv_to->logBF = get_logBF(cv_to->index,cv_to->size);
    //    cout << cv_to->logBF << endl;
    pFrom.AddToLogWeight(cv_to->logBF - prev_logBF);
  }
}

void fMove2(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng) // k-step look-ahead filter
{
  cv_state * cv_to = pFrom.GetValuePointer();

  /*
    double prev_weight_for_display = exp(pFrom.GetLogWeight());
    cout << setprecision(2) << "Particle " << particle_id << ": ";
    particle_id++;
  */
  
  if (cv_to->size < lTime - 1 || cv_to->size >= lIterates) {  
  }
  

  else {
    // build a subGBT that contains the variables in the current particle
    build_subGBT(*cv_to);    

    //  if (pRng->UniformS() < compute_rho(cv_to->size)) {  
    if (log(unifRand()) < models[0][1]) {
      pFrom.AddToLogWeight(log(compute_rho(cv_to->size,rho_method)) - models[0][1]);
    }
    
    else {
      // double u = pRng->UniformS();
      double u = unifRand();
      double prev_logBF;
      int i;
      //      double incr = 1.0/(p-cv_to->size);
      double cum_prob_curr = 0;
      
      double lambda_curr = 0;
      for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	lambda_curr = get_add_prob(cv_to,i);
	cum_prob_curr += lambda_curr;
      }
      


      cv_to->index = make_child_index(cv_to->index,i-1);
      cv_to->size += 1;
      prev_logBF = cv_to->logBF;
      
      /*
	cout << " stopping=" << exp(models[0][1]) << "; X" << i << " added, with prob = " << lambda_curr << " (";
	for (i=0; i< MAXVAR && cv_to->index.var[i]!=0; i++) {
	cout << cv_to->index.var[i] << ",";
	} 
	cout << ") ";
      */
      
      cv_to->logBF = get_logBF(cv_to->index,cv_to->size);
      //    cout << cv_to->logBF << endl;
      
      pFrom.AddToLogWeight(log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1])) - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF);
      
      //      cout << log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1])) - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF << endl;
      
      // cout << "rho_ratio=" << exp(log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1]))) << " " << "lambda_ratio=" << exp(- log(p-(cv_to->size-1))- log(lambda_curr)); 
    }

  }
  /*
    cout << "\t weight: " << prev_weight_for_display << "-->" << exp(pFrom.GetLogWeight());
    if (exp(pFrom.GetLogWeight())/prev_weight_for_display > 100) {
    cout << "***";
    }
    cout << endl;
  */
}


void fMove3(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng) // prior stopping + k-step selection
// This filter use the prior for proposing stopping and the k-step look-ahead to proposing selection
{
  cv_state * cv_to = pFrom.GetValuePointer();

  if (cv_to->size < lTime - 1 || cv_to->size >= lIterates ) {  
  }
  else {
    if (cv_to->size == MAXVAR || unifRand() < compute_rho(cv_to->size,rho_method)) {
    }
    
    else {
      // build a subGBT that contains the variables in the current particle
      build_subGBT(*cv_to);    
      
      // double u = pRng->UniformS();
      double u = unifRand();
      double prev_logBF;
      int i;
      //      double incr = 1.0/(p-cv_to->size);
      double cum_prob_curr = 0;
      
      double lambda_curr = 0;
      for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	lambda_curr = get_add_prob(cv_to,i);
	cum_prob_curr += lambda_curr;
      }
      
      // for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	
      // 	if (is_not_in_model(cv_to->index,i)) {
      // 	  cum_prob_curr += incr;
	  
      // 	}
      // }
      
      cv_to->index = make_child_index(cv_to->index,i-1);
      cv_to->size += 1;
      prev_logBF = cv_to->logBF;
      
      
      
      cv_to->logBF = get_logBF(cv_to->index,cv_to->size);
      //    cout << cv_to->logBF << endl;
      
      pFrom.AddToLogWeight( 0 - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF);

    }


  }
}

void fMove4(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng) // k-step stopping + prior selection
{
  cv_state * cv_to = pFrom.GetValuePointer();

  if (cv_to->size < lTime - 1 || cv_to->size >= lIterates) {  
  }
  
  //  else if (pRng->UniformS() < compute_rho(cv_to->size)) {
  else {
    // build a subGBT that contains the variables in the current particle
    build_subGBT(*cv_to);    

  
    if (log(unifRand()) < models[0][1]) {
      pFrom.AddToLogWeight(log(compute_rho(cv_to->size,rho_method)) - models[0][1]);
      //cout << "Stop! " << cv_to->size << "," << log(compute_rho(cv_to->size,rho_method)) << endl;
    }
    
    else {
      // double u = pRng->UniformS();
      double u = unifRand();
      double prev_logBF;
      int i;
      double incr = 1.0/(p-cv_to->size);
      double cum_prob_curr = 0;
      
      
      for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	
	if (is_not_in_model(cv_to->index,i)) {
	  cum_prob_curr += incr;
	  
	}
      }
      
      cv_to->index = make_child_index(cv_to->index,i-1);
      cv_to->size += 1;
      prev_logBF = cv_to->logBF;
      
      cv_to->logBF = get_logBF(cv_to->index,cv_to->size);
      //    cout << cv_to->logBF << endl;
      pFrom.AddToLogWeight(log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1])) + cv_to->logBF - prev_logBF);
    }
  }
}

void fMove5(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng) // pre-specified stopping + k-step selection
// This filter use the prespecified rho for proposing stopping and the k-step look-ahead to proposing selection
{
  cv_state * cv_to = pFrom.GetValuePointer();

  if (cv_to->size < lTime - 1 || cv_to->size >= lIterates ) {  
  }
  else {
    if (cv_to->size == MAXVAR || unifRand() < rho_prop) {
      pFrom.AddToLogWeight(log(compute_rho(cv_to->size,rho_method)) - log(rho_prop));
    }
    
    else {
      // build a subGBT that contains the variables in the current particle
      build_subGBT(*cv_to);    
      
      // double u = pRng->UniformS();
      double u = unifRand();
      double prev_logBF;
      int i;
      //      double incr = 1.0/(p-cv_to->size);
      double cum_prob_curr = 0;
      
      double lambda_curr = 0;
      for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	lambda_curr = get_add_prob(cv_to,i);
	cum_prob_curr += lambda_curr;
      }
      
      // for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	
      // 	if (is_not_in_model(cv_to->index,i)) {
      // 	  cum_prob_curr += incr;
	  
      // 	}
      // }
      
      cv_to->index = make_child_index(cv_to->index,i-1);
      cv_to->size += 1;
      prev_logBF = cv_to->logBF;
      
      
      
      cv_to->logBF = get_logBF(cv_to->index,cv_to->size);
      //    cout << cv_to->logBF << endl;
      
      pFrom.AddToLogWeight( log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-rho_prop) - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF);

    }

    
  }
}



void fMove6(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng) // SSS style proposal
// we give the equal selection probabilities to the top MAXVAR variables
{
  cv_state * cv_to = pFrom.GetValuePointer();

  
  if (cv_to->size < lTime - 1 || cv_to->size >= lIterates) {  
  }
  

  else {
    // build a subGBT that contains the variables in the current particle
    build_subGBT(*cv_to);    

    //  if (pRng->UniformS() < compute_rho(cv_to->size)) {  
    if (log(unifRand()) < models[0][1]) {
      pFrom.AddToLogWeight(log(compute_rho(cv_to->size,rho_method)) - models[0][1]);
    }
    
    else {
      // double u = pRng->UniformS();

      vector<ADD_PROB_TYPE> add_prob_vec(p);
      double u = unifRand();
      double prev_logBF;
      int i;
      int var_curr=0;
      double sum_top_add_probs;
      double cum_prob_curr = 0;
      double lambda_curr = 0;
      double top_incr;
      int n_top_vars_curr= max(1,min(n_top_vars,min(p,MAXVAR))-cv_to->size);
      
      for(i = 0; i < p; i++  ) {
	add_prob_vec[i] = make_pair(i,get_add_prob(cv_to,i));
      }
      
      sort(add_prob_vec.begin(),add_prob_vec.end(),comparator); // sort the add probabilities in desending order
      
      sum_top_add_probs=0;
      for(i = 0; i < n_top_vars_curr; i++) {
	sum_top_add_probs += add_prob_vec[i].second;
      }
      top_incr = sum_top_add_probs/n_top_vars_curr;
      // cout << top_incr << endl;

      cum_prob_curr = 0;
      for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	var_curr = add_prob_vec[i].first;
	
	
	if (i < n_top_vars_curr) lambda_curr = top_incr;
	else lambda_curr = add_prob_vec[i].second;
	

	// lambda_curr = add_prob_vec[i].second;

	cum_prob_curr += lambda_curr;
      }
      
      cv_to->index = make_child_index(cv_to->index,var_curr);
      cv_to->size += 1;
      prev_logBF = cv_to->logBF;
      
      /*
	cout << " stopping=" << exp(models[0][1]) << "; X" << i << " added, with prob = " << lambda_curr << " (";
	for (i=0; i< MAXVAR && cv_to->index.var[i]!=0; i++) {
	cout << cv_to->index.var[i] << ",";
	} 
	cout << ") ";
      */
      
      cv_to->logBF = get_logBF(cv_to->index,cv_to->size);
      //    cout << cv_to->logBF << endl;
      
      pFrom.AddToLogWeight(log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1])) - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF);
      
      //      cout << log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1])) - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF << endl;
      
      // cout << "rho_ratio=" << exp(log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1]))) << " " << "lambda_ratio=" << exp(- log(p-(cv_to->size-1))- log(lambda_curr)); 
    }

  }
  /*
    cout << "\t weight: " << prev_weight_for_display << "-->" << exp(pFrom.GetLogWeight());
    if (exp(pFrom.GetLogWeight())/prev_weight_for_display > 100) {
    cout << "***";
    }
    cout << endl;
  */
}

void fMove7(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng) // SSS style proposal
// we give the equal selection probabilities to the top MAXVAR variables
{
  cv_state * cv_to = pFrom.GetValuePointer();

  
  if (cv_to->size < lTime - 1 || cv_to->size >= lIterates) {  
  }
  

  else {
    // build a subGBT that contains the variables in the current particle
    build_subGBT(*cv_to);    

    //  if (pRng->UniformS() < compute_rho(cv_to->size)) {  
    if (log(unifRand()) < models[0][1]) {
      pFrom.AddToLogWeight(log(compute_rho(cv_to->size,rho_method)) - models[0][1]);
    }
    
    else {
      // double u = pRng->UniformS();

      vector<ADD_PROB_TYPE> add_prob_vec(p);
      double u = unifRand();
      double prev_logBF;
      int i;
      int var_curr=0;
      double sum_top_add_probs;
      double cum_prob_curr = 0;
      double lambda_curr = 0;
      double top_incr;
      // int n_top_vars_curr= min(n_top_vars,min(p,MAXVAR)-cv_to->size);
      
      for(i = 0; i < p; i++  ) {
	add_prob_vec[i] = make_pair(i,get_add_prob(cv_to,i));
      }
      
      sort(add_prob_vec.begin(),add_prob_vec.end(),comparator); // sort the add probabilities in descending order
      
      cum_prob_curr = 0;
      for(i = 0; i < p && cum_prob_curr < u; i++  ) {
	var_curr = add_prob_vec[i].first;
	lambda_curr = add_prob_vec[i].second;
	cum_prob_curr += lambda_curr;
      }
      
      cv_to->index = make_child_index(cv_to->index,var_curr);
      cv_to->size += 1;
      prev_logBF = cv_to->logBF;
      
      /*
	cout << " stopping=" << exp(models[0][1]) << "; X" << i << " added, with prob = " << lambda_curr << " (";
	for (i=0; i< MAXVAR && cv_to->index.var[i]!=0; i++) {
	cout << cv_to->index.var[i] << ",";
	} 
	cout << ") ";
      */
      
      cv_to->logBF = get_logBF(cv_to->index,cv_to->size);
      //    cout << cv_to->logBF << endl;
      
      pFrom.AddToLogWeight(log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1])) - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF);
      
      //      cout << log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1])) - log(p-(cv_to->size-1))- log(lambda_curr) + cv_to->logBF - prev_logBF << endl;
      
      // cout << "rho_ratio=" << exp(log(1-compute_rho(cv_to->size-1,rho_method)) - log(1-exp(models[0][1]))) << " " << "lambda_ratio=" << exp(- log(p-(cv_to->size-1))- log(lambda_curr)); 
    }

  }
  /*
    cout << "\t weight: " << prev_weight_for_display << "-->" << exp(pFrom.GetLogWeight());
    if (exp(pFrom.GetLogWeight())/prev_weight_for_display > 100) {
    cout << "***";
    }
    cout << endl;
  */
}






double get_add_prob(cv_state *model_ptr, int i) { //get the splitting probability of dimension i

  INDEX_TYPE I = model_ptr->index;
  int level = model_ptr->size;

  INDEX_TYPE I_sub = init_index(p-level,0);

  double *node = get_node(I_sub,0);

  if (is_not_in_model(I,i)) { // if the ith predictor is not in the model

    int i_sub = get_i_sub(I,i);
    double loglambda0 = (-1.0) * log(p - level);
    return (1.0 - compute_rho(level,rho_method)) * exp( loglambda0 + get_child(I_sub,i_sub,0)[0] - node[0] - log(1 -exp(node[1])) ); 
  } else {
    return 0.0;
  } // if ith predictor is already in the model then return 0;
}

inline int get_i_sub(INDEX_TYPE I, int i) {
  int i_sub = i;
  for (int j = 0; I.var[j] > 0 && I.var[j] < i + 1 && j < MAXVAR; j++) {
    i_sub--;
  }
  return i_sub;
}


void build_subGBT(const cv_state & root_model) { // return a vector of length p+1 
  // the first p elements being the posterior selection probabilities
  // the last element is the posterior logrho

  int i, level;
  double dlogrho0;
  double dlogrho1;
  double d_logphi;
  double loglambda0;
  double m_stop,m_divide;
  int p_sub = p - root_model.size;  
  int k_sub = min(kstep, min(p_sub,MAXVAR-root_model.size));
  INDEX_TYPE I_global,I_sub;
  ulong count;
  double *NODE_CURR;
  //  models = new double*[k_sub+1];
  // modelscount = new unsigned int[k_sub+1];

    
  for (level=k_sub; level >=0; level--) {

    count = 0;
    // I = init_index_subGBT(root_model.index,level);
    I_sub = init_index(p_sub,level);
    loglambda0 = (-1.0) * log(p_sub - level);

    dlogrho0 = log(compute_rho(root_model.size + level, rho_method));
    dlogrho1 = log(1-exp(dlogrho0));

    // while (count < modelscount[level]) { 
    while (count < Choose(p_sub,level)) { 
      I_global = make_global_index(root_model.index,I_sub,level);

      //NODE_CURR = get_node_subGBT(root_model.index,I,level);
      NODE_CURR = get_node(I_sub,level);

      NODE_CURR[0] = get_logBF(I_global,root_model.size+level); 
      // at this point NODE_CURR[0] contains the logBF for model I_global
      // In the end, NODE_CURR[0] will contain phi
      
      if (level == k_sub) {
	// do nothing so NODE_CURR[0] is equal to \phi
	NODE_CURR[1] = 0; // NODE_CURR[1] = posterior logrho
      }
      
      else {
	m_divide = -DBL_MAX;
	for (i=0; i < p_sub; i++) { // for each dimension i
	  if (is_not_in_model(I_sub,i)) { // if the ith predictor is addable

	    d_logphi = get_child(I_sub,i,level)[0];

	    if (m_divide == -DBL_MAX) {
	      m_divide = d_logphi;
	    } else {
	      m_divide = log_exp_x_plus_exp_y (m_divide, d_logphi);
	    }    
	  }	  
	}
	

	m_divide += dlogrho1 + loglambda0; 
	m_stop = dlogrho0 + NODE_CURR[0];
	NODE_CURR[0] = log_exp_x_plus_exp_y(m_stop, m_divide);
	NODE_CURR[1] = m_stop - NODE_CURR[0];
      }     
      
      I_sub = get_next_node(I_sub,p_sub,level); count++;
    }
  }  
}

INDEX_TYPE make_global_index(INDEX_TYPE I_root, INDEX_TYPE I_sub, int level) {
  
  int i,j=0;
  INDEX_TYPE I_global = I_root;

  for (i =0; I_sub.var[i] > 0 && i < level; i++) {
    for (; I_root.var[j] - j - 1 < I_sub.var[i] && I_root.var[j] > 0 && j < MAXVAR; j++);
    // after this loop j is equal to the number of variables before I_sub.var[i]

    if (j == MAXVAR){
      cout << "Warning: Variable cannot be added, reached maximum model size allowed! This should not happen!" << endl;
    }
    else {
      I_global = make_child_index(I_global, I_sub.var[i]+j-1); // Add variable I_sub.var[i] + j
    }
  }

  return I_global;
}


inline INDEX_TYPE make_child_index(INDEX_TYPE& I, unsigned short add_dim) {
    INDEX_TYPE child_index = I;
    unsigned short data = add_dim+1; 
    int i = 0;
    while (i<MAXVAR) {
        while (child_index.var[i] >0 && data > child_index.var[i] ) {
            i++;
        }
        if (child_index.var[i] == 0) {
            child_index.var[i] = data;
            i = MAXVAR;
        } else {
            unsigned short swap = child_index.var[i];
            child_index.var[i] = data; data= swap;
            i++;
        }
    }
    return child_index;
}


double *get_child(INDEX_TYPE& I, int i,int level) {
     INDEX_TYPE child_index = make_child_index(I,i);
     return &models[level+1][(get_node_index(child_index,level+1))];
}

inline INDEX_TYPE get_next_node(INDEX_TYPE& I, int n, int k) {
  //print_index(I,k); cout << "--->";
    INDEX_TYPE node = I;
    int i = k-1; int j = n-1;
    while (i>=0 && node.var[i] == j+1) {i--;j--;}
    if (i < 0) { //reach the end of nodes
        node.index=0; //invalid node
    } else {
        node.var[i] += 1;
        for (j=i+1;j<k;j++) {
            node.var[j] = node.var[i]+ j-i;
        }
    }

    //print_index(node,k); cout << endl;
    return node;
    
}

double *get_node(INDEX_TYPE& I, int level) {
    return &models[level][(get_node_index(I,level))];
}

unsigned int get_node_index(INDEX_TYPE& I,int level) {
     unsigned long long  r = 0;
 
     unsigned long long numerator = 1;
     unsigned long long denominator = 1;
     
     for (int i = 0; i < level; i++) {
       numerator = 1;
       denominator *= (i+1);     
       for (int j = 1; j <= i+1; j++) {
     	 numerator *= I.var[i]-j;
      }
       r += numerator / denominator;
     }

     return r*NUMNODEVAR;
}

long fSelect(long lTime, const smc::particle<cv_state> & part,smc::rng *pRng)
{

  // return 0; // bootstrap filter 
  // return 1; // k-step look-ahead filter
  return (unifRand() > p_mix); // a half-half mixture of the two filters. 

}
