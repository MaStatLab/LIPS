#ifdef RLIB
    #include "hyp2f1.c"
#endif
double unifRand() {
  return rand() / ( (double) RAND_MAX + 1);
}
double ddot(double *a, double *b, int p) {
    double r = 0.0;
    for (int i = 0; i < p;i++) {
        r+=a[i]*b[i];
    }
    return r;
}

#ifdef MATLAB
    void cholreg(double *XtX, double *XtY,  double *coefficients, double *mse,  int p, int n)
    {
        // On entry *coefficients undefined, which is over written with the OLS estimates 
        // On entry MSE = Y'Y 
        char Uplo = 'U';
        ptrdiff_t info,col=1;
        ptrdiff_t pp = p;
        double *A2 = new double[p*p];
        memcpy(A2, XtX, p*p*sizeof(double));
        memcpy(coefficients, XtY, p*sizeof(double));
        dposv(&Uplo, &pp,&col,A2,&pp,coefficients,&pp,&info);
        double ete = ddot(XtY,coefficients,p);
        *mse = (*mse - ete)/((double) (n - p));
        delete [] A2;
    } 
#else
    //In-place Chol Descoposition
    //no extra space are used and both upper and lower 
    //matrices are saved in the same return matrix
    int Chol_Inplace(double y[], int n) {
        double temp;
        int i,j,k;
        for (i=0; i<n;i++) {
            for (j = 0; j < i;j++) {
                temp = y[i+j*n];
                for (k=0; k <j;k++) {
                    temp-=y[i+k*n] * y[j+k*n];
                }
                y[i+j*n] = temp / y[j+j*n];
                y[j+i*n] = y[i+j*n];
            }
            temp = y[i+i*n];
            if (temp == 0.0) return -1;
            for (k=0; k <i;k++) {
                temp-=y[i+k*n] * y[i+k*n];
            }
            y[i+i*n] = sqrt(temp);
        }
        return 0;
    }

    int Lower_Triangular_Solve(double *L, double B[], double x[], int n) {
        int i, k;

        //         Solve the linear equation Lx = B for x, where L is a lower
        //         triangular matrix.                                      

        for (k = 0; k < n; L += n, k++) {
            if (*(L + k) == 0.0) return -1;           // The matrix L is singular
            x[k] = B[k];
            for (i = 0; i < k; i++) x[k] -= x[i] * *(L + i);
            x[k] /= *(L + k);
        }
        return 0;
    }
    int Upper_Triangular_Solve(double *U, double B[], double x[], int n) {
        int i, k;

        //         Solve the linear equation Ux = B for x, where U is an upper
        //         triangular matrix.                                      

        for (k = n-1, U += n * (n - 1); k >= 0; U -= n, k--) {
            if (*(U + k) == 0.0) return -1;           // The matrix U is singular
            x[k] = B[k];
            for (i = k + 1; i < n; i++) x[k] -= x[i] * *(U + i);
            x[k] /= *(U + k);
        }
        return 0;
    }

    void cholreg(double *XtX, double *XtY,  double *coefficients, double *mse,  int p, int n)
    {
        // On entry *coefficients undefined, which is over written with the OLS estimates 
        // On entry MSE = Y'Y 
        double *work = new double[p*p];
        double *y0 = new double[p];
        memcpy(work,XtX,p*p*sizeof(double));
        Chol_Inplace(work,p);
        Lower_Triangular_Solve(work, XtY, y0, p);
        Upper_Triangular_Solve(work, y0, coefficients, p);

        double ete = ddot(XtY,coefficients,p);
        *mse = (*mse - ete)/((double) (n - p));
        free(work);
        free(y0);

    } 
#endif
double log_exp_x_plus_exp_y(double x, double y) {

    double result;
    if ( x - y >= 100 ) result = x;
    else if ( x - y <= -100 ) result = y;
    else {
      if (x > y) {
      result = y + log( 1 + exp(x-y) );
      }
      else result = x + log( 1 + exp(y-x) );
    }
    return result;
}

double logBF_hyperGprior_laplace(double R2, int n, int p, double alpha)
 {  
   /* R2 = usual coefficient of determination
      n = sample size
      p = number of rank of X (including intercept)
      alpha = prior hyperparameter
      n and p are adjusted by subtrating off one 
      because of the flat prior on the intercept
   */

   double lognc, ghat, dn, dp, logmarg,sigmahat;
   
    dn = (double) n - 1.0;
    dp = (double) p - 1.0;
/*  Laplace approximation in terms of exp(tau) = g  */
/*  Agrees with Mathematica but not paper  */
    ghat = (-4.+ alpha + dp + (2. - dn)*R2 - 
	    sqrt(-8.*(-2. + alpha + dp)*(-1.0 + R2) + (-4. + alpha + dp + (2.-dn)* R2)*(-4. + alpha + dp + (2.-dn)* R2)))/(2.*(-2. + alpha + dp)*(-1. + R2)); 

    //if (ghat <= 0.0)  { Rprintf("ERROR: In Laplace approximation to  logmarg,  ghat =  %f  R2 = %f p = %d  n = %d\n", ghat, R2, p,n);}
  
   
    /*  Substitute ghat (unconstrained) into sigma, simplify, then plug in ghat
	Seems to perform better */

        
    sigmahat =1.0/(-ghat*(dn - alpha - dp)/(2.*(1. + ghat)*(1.+ghat)) +
                    dn*(ghat*(1. - R2))/(2.*(1.+ghat*(1.-R2))*(1.+ghat*(1.-R2)))); 

    //if (sigmahat <= 0 ) Rprintf("ERROR in LAPLACE APPROXIMATION to logmarg sigmhat = %f, ghat =  %f  R2 = %f p = %d  n = %d\n", sigmahat, ghat, R2, p,n); 
    lognc = log(alpha/2.0 - 1.0);
    logmarg = lognc + 
              .5*( log(2.0*PI) 
                     - (dp + alpha)*log(1.0 + ghat)
	             -  dn*log(1.0-(ghat/(1.0 + ghat))*R2)
	             + log(sigmahat)) + log(ghat);
  if (p == 1) logmarg = 0.0;
  return(logmarg);
}

double logBF_gprior( double Rsquare, int n,  int p, double g)
  {  double logmargy;
  /* log Marginal under reference prior for phi, intercept and
     g prior for regression coefficients 
     log marginal for null model is 0 */
    logmargy =  .5*(log(1.0 + g) * ((double) (n - p))  - log(1.0 + g * (1.0- Rsquare)) * ((double)(n - 1)));
  if (p == 1) logmargy = 0.0;
  return(logmargy);
  }

#ifdef RLIB
extern double hyp2f1(double a, double b, double c, double x);

double logBF_hyperGprior(double R2, int n, int p, double alpha)
{  double logmargy,  a, b, c, z1, hf1;

  a = (double)  (n - 1) /2.0;
  b = 1.0;
  c = ((double) p - 1.0  + alpha )/2.0;
  z1 = R2;
  
  hf1 = hyp2f1(a, b, c, z1);
  if (p == 1) logmargy = 0.0;
  else logmargy = log(hf1) 
	 - log( (double) p - 1.0 + alpha - 2.0) + log(2.0) 
	 + log(alpha/2.0 - 1.0);
  if (! R_FINITE(logmargy))
    logmargy = logBF_hyperGprior_laplace(R2, n, p, alpha);
  return(logmargy);
}
#endif

