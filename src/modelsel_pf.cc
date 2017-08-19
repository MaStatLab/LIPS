#include <R.h>
#include <Rinternals.h>

// #include <Rcpp.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <list>
#include <map>
#define RLIB

// #include "modeltree_node.cpp"
// #include "gbt.cpp"


#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "smctc/include/smctc.hh"
#include "modelsel_pf_funcs.cc"


double yty = 0;
double SSY = 0;
int nobs = 0;
int p = 0;
int method = 0;
int lIterates;
int kstep;
double alpha = 0;
double *XtX;
double *XtY;
double *XtXwork;
double *XtYwork;
double *coefficients;
double rho0;
int rho_method;
double **models;
double p_mix;
double rho_prop;
int n_top_vars;

int particle_id;

extern "C" {

  SEXP modelsel_pf_C(SEXP Y, SEXP X, SEXP n_pred, SEXP T, SEXP k_step, SEXP Ralpha, SEXP Rcoef_prior, SEXP Rn_particle,SEXP Rrho0, SEXP Rrho_method, SEXP Rrho_prop, SEXP Rp_mix, SEXP Rresample_param,SEXP Rn_top_vars);
  SEXP modelsel_pf_hd_C(SEXP Y, SEXP X, SEXP n_pred, SEXP T, SEXP k_step, SEXP Ralpha, SEXP Rcoef_prior, SEXP Rn_particle,SEXP Rrho0, SEXP Rrho_method, SEXP Rrho_prop, SEXP Rp_mix, SEXP Rresample_param,SEXP Rn_top_vars);
  SEXP modelsel_pf_hd_pred_C(SEXP Y, SEXP X, SEXP n_pred, SEXP T, SEXP k_step, SEXP Ralpha, SEXP Rcoef_prior, SEXP Rn_particle, SEXP Rrho0, SEXP Rrho_method, SEXP Rrho_prop, SEXP Rp_mix, SEXP Rresample_param, SEXP Rn_top_vars);
};

SEXP modelsel_pf_C(SEXP Y, SEXP X, SEXP n_pred, SEXP T, SEXP k_step, SEXP Ralpha, SEXP Rcoef_prior, SEXP Rn_particle, SEXP Rrho0, SEXP Rrho_method, SEXP Rrho_prop, SEXP Rp_mix, SEXP Rresample_param, SEXP Rn_top_vars) {

  SEXP RXwork = PROTECT(duplicate(X)), RYwork = PROTECT(duplicate(Y));
  int nProtected = 2;

  double *Xwork, *Ywork;
  int n_part,lTime;
  int i,j,l;
  cv_state model_curr;

  double ybar;

  nobs = LENGTH(Y);
  p = INTEGER(n_pred)[0];
  lIterates = min(min(INTEGER(T)[0],p),MAXVAR);
  kstep = min(INTEGER(k_step)[0],lIterates);
  n_part = INTEGER(Rn_particle)[0];
  method = INTEGER(Rcoef_prior)[0];
  alpha = REAL(Ralpha)[0];
  rho0 = REAL(Rrho0)[0];
  rho_method = INTEGER(Rrho_method)[0];
  rho_prop = REAL(Rrho_prop)[0];
  p_mix = REAL(Rp_mix)[0];
  n_top_vars = min(min(MAXVAR,p),INTEGER(Rn_top_vars)[0]);
  double resample_param = REAL(Rresample_param)[0];

  // First we read in the data
  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);

  // these are used to store the entire data set
  XtX  = new double[p*p];
  XtY = new double[p];


  memset(XtX,0,p*p*sizeof(double));
  ybar = 0.0; SSY = 0.0; yty = 0.0;
  for (i=0; i <p; i++) { //since we do it only once, do this for now
    for (j=i; j < p; j++) {
      double sum = 0.0;
      for (l=0;l<nobs;l++) {
	sum += Xwork[i*nobs+l] * Xwork[j*nobs+l];
      }
      XtX[i*p+j] = sum; XtX[j*p+i] = sum;
    }
  }
  for (i =0; i < nobs; i++) {
    yty += Ywork[i]*Ywork[i];
    ybar += Ywork[i];
  }
  ybar = ybar/ (double) nobs;
  SSY = yty - (double) nobs* ybar *ybar;
  for (i=0; i <p; i++) {
      double sum = 0.0;
      for (j=0;j<nobs;j++) {
	sum += Xwork[i*nobs+j] * Ywork[j];
      }
      XtY[i] = sum;
  }
  // Up to this point the data are stored in XtX and XtY

  int temp = min(MAXVAR,lIterates);

  if (p_mix != 1) {

    models = new double*[kstep+1];

    for (i=0; i <= kstep; i++) {
      models[i] = new double[Choose(p,i)*NUMNODEVAR];
    }
  }

  try {
    // Run particle filter


    smc::sampler<cv_state> Sampler(n_part, SMC_HISTORY_NONE);
    if (p_mix == 1) {
      smc::moveset<cv_state> Moveset(fInitialise, fMove, NULL);
      Sampler.SetMoveSet(Moveset);
    }
    else if (p_mix == 0) {
      smc::moveset<cv_state> Moveset(fInitialise, fMove2, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 2) { // hybrid filter prior stopping k-step selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove3, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 3) { // hybrid filter k-step stopping prior selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove4, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 4) { // hybrid filter k-step stopping prior selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove5, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 5) { // SSS style proposal
      smc::moveset<cv_state> Moveset(fInitialise, fMove6, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 6) { // same as fMove2 but with sorted add probability vector
      smc::moveset<cv_state> Moveset(fInitialise, fMove7, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (0 < p_mix && p_mix < 1) {
      void (*pfMoves[])(long, smc::particle<cv_state> &,smc::rng*) = {fMove, fMove2};
      smc::moveset<cv_state> Moveset(fInitialise, fSelect,sizeof(pfMoves)/sizeof(pfMoves[0]),pfMoves,NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else {
      cerr << "Error: p.mix must be in [0,1], 2 or 3" << endl;
      exit(1);
    }

    Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, resample_param);
    Sampler.Initialise();


    for(lTime=1 ; lTime <= lIterates ; lTime++) {
      particle_id = 0;
      Sampler.Iterate();
      cout << "Step " << lTime << ": ESS = " << Sampler.GetESS() << endl;
    }


    // Return the results
    SEXP particle_models;
    SEXP particle_log_weights;

    PROTECT(particle_models = allocMatrix(INTSXP,n_part,p)); ++nProtected;
    PROTECT(particle_log_weights = allocVector(REALSXP,n_part)); ++nProtected;

    for(i = 0; i < n_part; i++) {

      model_curr = Sampler.GetParticleValue(i);

      for (j = 0; j < p; j++) {

	      INTEGER(particle_models)[i + n_part*j] = is_in_model(model_curr.index,j);

      }

      REAL(particle_log_weights)[i] = Sampler.GetParticleLogWeight(i);

    }

    SEXP ans;
    PROTECT(ans = allocVector(VECSXP,2)); ++nProtected;
    SET_VECTOR_ELT(ans,0,particle_models);
    SET_VECTOR_ELT(ans,1,particle_log_weights);



    delete [] XtX;
    delete [] XtY;

    if (p_mix != 1) {

      for (i=0; i <= kstep; i++) {
	delete [] models[i];
      }

      delete [] models;
    }

    UNPROTECT(nProtected);
    return ans;

  }

  catch(smc::exception e) {
    cerr << e;
    exit(e.lCode);
  }
}

SEXP modelsel_pf_hd_C(SEXP Y, SEXP X, SEXP n_pred, SEXP T, SEXP k_step, SEXP Ralpha, SEXP Rcoef_prior, SEXP Rn_particle, SEXP Rrho0, SEXP Rrho_method, SEXP Rrho_prop, SEXP Rp_mix, SEXP Rresample_param, SEXP Rn_top_vars) {

  SEXP RXwork = PROTECT(duplicate(X)), RYwork = PROTECT(duplicate(Y));
  int nProtected = 2;

  double *Xwork, *Ywork;
  int n_part,lTime;
  int i,j,l;
  cv_state model_curr;

  double ybar;

  nobs = LENGTH(Y);
  p = INTEGER(n_pred)[0];
  lIterates = min(min(INTEGER(T)[0],p),MAXVAR);
  kstep = min(INTEGER(k_step)[0],lIterates);
  n_part = INTEGER(Rn_particle)[0];
  method = INTEGER(Rcoef_prior)[0];
  alpha = REAL(Ralpha)[0];
  rho0 = REAL(Rrho0)[0];
  rho_method = INTEGER(Rrho_method)[0];
  rho_prop = REAL(Rrho_prop)[0];
  p_mix = REAL(Rp_mix)[0];
  n_top_vars = min(min(MAXVAR,p),INTEGER(Rn_top_vars)[0]);
  double resample_param = REAL(Rresample_param)[0];

  // First we read in the data
  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);

  // these are used to store the entire data set
  XtX  = new double[p*p];
  XtY = new double[p];

  memset(XtX,0,p*p*sizeof(double));
  ybar = 0.0; SSY = 0.0; yty = 0.0;
  for (i=0; i <p; i++) { //since we do it only once, do this for now
    for (j=i; j < p; j++) {
      double sum = 0.0;
      for (l=0;l<nobs;l++) {
	sum += Xwork[i*nobs+l] * Xwork[j*nobs+l];
      }
      XtX[i*p+j] = sum; XtX[j*p+i] = sum;
    }
  }
  for (i =0; i < nobs; i++) {
    yty += Ywork[i]*Ywork[i];
    ybar += Ywork[i];
  }
  ybar = ybar/ (double) nobs;
  SSY = yty - (double) nobs* ybar *ybar;
  for (i=0; i <p; i++) {
      double sum = 0.0;
      for (j=0;j<nobs;j++) {
	sum += Xwork[i*nobs+j] * Ywork[j];
      }
      XtY[i] = sum;
  }
  // Up to this point the data are stored in XtX and XtY

  int temp = min(MAXVAR,lIterates);

  if (p_mix != 1) {

    models = new double*[kstep+1];

    for (i=0; i <= kstep; i++) {
      models[i] = new double[Choose(p,i)*NUMNODEVAR];
    }
  }

  try {
    // Run particle filter


    smc::sampler<cv_state> Sampler(n_part, SMC_HISTORY_NONE);
    if (p_mix == 1) {
      smc::moveset<cv_state> Moveset(fInitialise, fMove, NULL);
      Sampler.SetMoveSet(Moveset);
    }
    else if (p_mix == 0) {
      smc::moveset<cv_state> Moveset(fInitialise, fMove2, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 2) { // hybrid filter prior stopping k-step selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove3, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 3) { // hybrid filter k-step stopping prior selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove4, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 4) { // hybrid filter k-step stopping prior selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove5, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 5) { // SSS style proposal
      smc::moveset<cv_state> Moveset(fInitialise, fMove6, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 6) { // same as fMove2 but with sorted add probability vector
      smc::moveset<cv_state> Moveset(fInitialise, fMove7, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (0 < p_mix && p_mix < 1) {
      void (*pfMoves[])(long, smc::particle<cv_state> &,smc::rng*) = {fMove, fMove2};
      smc::moveset<cv_state> Moveset(fInitialise, fSelect,sizeof(pfMoves)/sizeof(pfMoves[0]),pfMoves,NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else {
      cerr << "Error: p.mix must be in [0,1], 2 or 3" << endl;
      exit(1);
    }

    Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, resample_param);
    Sampler.Initialise();

    cout << "lIterates = " << lIterates << endl;

    for(lTime=1 ; lTime <= lIterates ; lTime++) {
      particle_id = 0;
      Sampler.Iterate();
      cout << "Step " << lTime << ": ESS = " << Sampler.GetESS() << endl;
    }

    SEXP particle_log_weights;
    SEXP avg_incl_probs;
    SEXP model_size;
    SEXP model_logBF;

    PROTECT(particle_log_weights = allocVector(REALSXP,n_part)); ++nProtected;
    PROTECT(avg_incl_probs = allocVector(REALSXP,p)); ++nProtected;
    PROTECT(model_size = allocVector(INTSXP,n_part)); ++nProtected;
    PROTECT(model_logBF = allocVector(REALSXP,n_part)); ++nProtected;

    for (j = 0; j < p; j++) REAL(avg_incl_probs)[j] = 0;

    double total_weight = 0;

    for(i = 0; i < n_part; i++) {

      total_weight += Sampler.GetParticleWeight(i);


      model_curr = Sampler.GetParticleValue(i);
      INTEGER(model_size)[i] = 0;

      for (j = 0; j < min(MAXVAR,lIterates); j++) {
	      if (model_curr.index.var[j] >0) {
	        REAL(avg_incl_probs)[model_curr.index.var[j]-1] += Sampler.GetParticleWeight(i);
	        INTEGER(model_size)[i] += 1;
	      }
      }

      REAL(particle_log_weights)[i] = Sampler.GetParticleLogWeight(i);
      REAL(model_logBF)[i] = model_curr.logBF;
    }

    for (j = 0; j < p; j++) {
      REAL(avg_incl_probs)[j] /= total_weight;
    }

    SEXP ans;
    PROTECT(ans = allocVector(VECSXP,4)); ++nProtected;
    SET_VECTOR_ELT(ans,0,particle_log_weights);
    SET_VECTOR_ELT(ans,1,avg_incl_probs);
    SET_VECTOR_ELT(ans,2,model_size);
    SET_VECTOR_ELT(ans,3,model_logBF);

    delete [] XtX;
    delete [] XtY;

    if (p_mix != 1) {

      for (i=0; i <= kstep; i++) {
    	delete [] models[i];
      }

      delete [] models;
    }

    UNPROTECT(nProtected);

    return ans;

  }

  catch(smc::exception e) {
    cerr << e;
    exit(e.lCode);
  }
}




SEXP modelsel_pf_hd_pred_C(SEXP Y, SEXP X, SEXP n_pred, SEXP T, SEXP k_step, SEXP Ralpha, SEXP Rcoef_prior, SEXP Rn_particle, SEXP Rrho0, SEXP Rrho_method, SEXP Rrho_prop, SEXP Rp_mix, SEXP Rresample_param, SEXP Rn_top_vars) {

  SEXP RXwork = PROTECT(duplicate(X)), RYwork = PROTECT(duplicate(Y));
  int nProtected = 2;

  double *Xwork, *Ywork;
  int n_part,lTime;
  int i,j,l;
  cv_state model_curr;

  double ybar;

  nobs = LENGTH(Y);
  p = INTEGER(n_pred)[0];
  lIterates = min(min(INTEGER(T)[0],p),MAXVAR);
  kstep = min(INTEGER(k_step)[0],lIterates);
  n_part = INTEGER(Rn_particle)[0];
  method = INTEGER(Rcoef_prior)[0];
  alpha = REAL(Ralpha)[0];
  rho0 = REAL(Rrho0)[0];
  rho_method = INTEGER(Rrho_method)[0];
  rho_prop = REAL(Rrho_prop)[0];
  p_mix = REAL(Rp_mix)[0];
  n_top_vars = min(min(MAXVAR,p),INTEGER(Rn_top_vars)[0]);
  double resample_param = REAL(Rresample_param)[0];

  // First we read in the data
  Ywork = REAL(RYwork);
  Xwork = REAL(RXwork);

  // these are used to store the entire data set
  XtX  = new double[p*p];
  XtY = new double[p];

  memset(XtX,0,p*p*sizeof(double));
  ybar = 0.0; SSY = 0.0; yty = 0.0;
  for (i=0; i <p; i++) { //since we do it only once, do this for now
    for (j=i; j < p; j++) {
      double sum = 0.0;
      for (l=0;l<nobs;l++) {
	sum += Xwork[i*nobs+l] * Xwork[j*nobs+l];
      }
      XtX[i*p+j] = sum; XtX[j*p+i] = sum;
    }
  }
  for (i =0; i < nobs; i++) {
    yty += Ywork[i]*Ywork[i];
    ybar += Ywork[i];
  }
  ybar = ybar/ (double) nobs;
  SSY = yty - (double) nobs* ybar *ybar;
  for (i=0; i <p; i++) {
      double sum = 0.0;
      for (j=0;j<nobs;j++) {
	sum += Xwork[i*nobs+j] * Ywork[j];
      }
      XtY[i] = sum;
  }
  // Up to this point the data are stored in XtX and XtY

  int temp = min(MAXVAR,lIterates);

  if (p_mix != 1) {

    models = new double*[kstep+1];

    for (i=0; i <= kstep; i++) {
      models[i] = new double[Choose(p,i)*NUMNODEVAR];
    }
  }

  try {
    // Run particle filter


    smc::sampler<cv_state> Sampler(n_part, SMC_HISTORY_NONE);
    if (p_mix == 1) {
      smc::moveset<cv_state> Moveset(fInitialise, fMove, NULL);
      Sampler.SetMoveSet(Moveset);
    }
    else if (p_mix == 0) {
      smc::moveset<cv_state> Moveset(fInitialise, fMove2, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 2) { // hybrid filter prior stopping k-step selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove3, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 3) { // hybrid filter k-step stopping prior selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove4, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 4) { // hybrid filter k-step stopping prior selection
      smc::moveset<cv_state> Moveset(fInitialise, fMove5, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 5) { // SSS style proposal
      smc::moveset<cv_state> Moveset(fInitialise, fMove6, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (p_mix == 6) { // same as fMove2 but with sorted add probability vector
      smc::moveset<cv_state> Moveset(fInitialise, fMove7, NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else if (0 < p_mix && p_mix < 1) {
      void (*pfMoves[])(long, smc::particle<cv_state> &,smc::rng*) = {fMove, fMove2};
      smc::moveset<cv_state> Moveset(fInitialise, fSelect,sizeof(pfMoves)/sizeof(pfMoves[0]),pfMoves,NULL);
      Sampler.SetMoveSet(Moveset);
    }

    else {
      cerr << "Error: p.mix must be in [0,1], 2 or 3" << endl;
      exit(1);
    }

    Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, resample_param);
    Sampler.Initialise();

    cout << "lIterates = " << lIterates << endl;

    for(lTime=1 ; lTime <= lIterates ; lTime++) {
      particle_id = 0;
      Sampler.Iterate();
      cout << "Step " << lTime << ": ESS = " << Sampler.GetESS() << endl;
    }

    // Return the results
    SEXP particle_models;
    SEXP particle_log_weights;
    SEXP avg_incl_probs;
    SEXP model_size;
    SEXP model_logBF;

    PROTECT(particle_models = allocMatrix(INTSXP,n_part,min(MAXVAR,lIterates))); ++nProtected;

    PROTECT(particle_log_weights = allocVector(REALSXP,n_part)); ++nProtected;
    PROTECT(avg_incl_probs = allocVector(REALSXP,p)); ++nProtected;
    PROTECT(model_size = allocVector(INTSXP,n_part)); ++nProtected;
    PROTECT(model_logBF = allocVector(REALSXP,n_part)); ++nProtected;

    for (j = 0; j < p; j++) REAL(avg_incl_probs)[j] = 0;

    double total_weight = 0;

    for(i = 0; i < n_part; i++) {

      total_weight += Sampler.GetParticleWeight(i);


      model_curr = Sampler.GetParticleValue(i);
      INTEGER(model_size)[i] = 0;

      for (j = 0; j < min(MAXVAR,lIterates); j++) {
	      if (model_curr.index.var[j] >0) {

	        REAL(avg_incl_probs)[model_curr.index.var[j]-1] += Sampler.GetParticleWeight(i);
	        INTEGER(model_size)[i] += 1;
	      }
	      INTEGER(particle_models)[i + n_part*j] = model_curr.index.var[j];

      }

      REAL(particle_log_weights)[i] = Sampler.GetParticleLogWeight(i);
      REAL(model_logBF)[i] = model_curr.logBF;

    }

    for (j = 0; j < p; j++) {
      REAL(avg_incl_probs)[j] /= total_weight;
    }

    SEXP ans;
    PROTECT(ans = allocVector(VECSXP,5)); ++nProtected;
    SET_VECTOR_ELT(ans,0,particle_log_weights);
    SET_VECTOR_ELT(ans,1,avg_incl_probs);
    SET_VECTOR_ELT(ans,2,model_size);
    SET_VECTOR_ELT(ans,3,model_logBF);
    SET_VECTOR_ELT(ans,4,particle_models);

    delete [] XtX;
    delete [] XtY;

    if (p_mix != 1) {

      for (i=0; i <= kstep; i++) {
    	delete [] models[i];
      }

      delete [] models;
    }

    UNPROTECT(nProtected);

    return ans;

  }

  catch(smc::exception e) {
    cerr << e;
    exit(e.lCode);
  }
}
