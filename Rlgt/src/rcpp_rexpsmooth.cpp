#include <Rcpp.h>
using namespace Rcpp;

int grid_sample(double u, NumericVector Lprop, double Lprop_min)
{
  int nProps = Lprop.size();
  int i;
  
  // Renormalize to probabilities
  double pSum = 0;
  for (i = 0; i<nProps; i++)
  {
    Lprop[i] = exp(-(Lprop[i] - Lprop_min));
    pSum = pSum + Lprop[i];
  }
  double us = u*pSum;
  
  // Now sample ...
  double acc = 0;
  for (i = 0; i<nProps; i++)
  {
    if (us > acc && us <= acc+Lprop[i])
    {
      return i;
    }
    acc = acc + Lprop[i];
  }
  
  return i;
}

// [[Rcpp::export]]
NumericVector rcpp_GridSampleTauPhi(NumericMatrix T, NumericVector u, NumericVector chi2, NumericVector e, NumericVector logl, NumericVector omega2, NumericVector nu)
{
  //
  int i, j;
  int nProps = T.rows();
  int n      = e.size();
  double phi, tau;
  double v2;
  double Lprop_min = 1e100;
  
  NumericVector Lprop(nProps);

  // Compute negative log-posteriors for each proposal
  for (i = 0; i<nProps; i++)
  {
    Lprop[i] = 0;
    // Get proposed parameter values
    phi = T(i,0);
    tau = T(i,1);

    for (j = 0; j<n; j++)
    {
      v2 = chi2[0]*(phi*phi + (1-phi)*(1-phi)*exp(logl[j]*2.0*tau));
      Lprop[i] = Lprop[i] + (nu[0]+1.0)*0.5*log(1 + e[j]*e[j]/v2/nu[0]) + 0.5*log(v2);
      //Lprop[i] = Lprop[i] + 0.5*e[j]*e[j]/v2 + log(v2)*0.5;
    }
    if (Lprop[i] < Lprop_min)
    {
      Lprop_min = Lprop[i];
    }
  }
  
  // // Renormalize to probabilities
  // double pSum = 0;
  // for (i = 0; i<nProps; i++)
  // {
  //   Lprop[i] = exp(-(Lprop[i] - Lprop_min));
  //   pSum = pSum + Lprop[i];
  // }
  // double us = u[0]*pSum;
  // 
  // // Now sample ...
  // double acc = 0;
  // for (i = 0; i<nProps; i++)
  // {
  //   if (us > acc && us <= acc+Lprop[i])
  //   {
  //     break;
  //   }
  //   acc = acc + Lprop[i];
  // }
  
  //i = grid_sample(u[0], Lprop, Lprop_min)+1;
  NumericVector theta = T.row(grid_sample(u[0], Lprop, Lprop_min));
  return theta;
}

// [[Rcpp::export]]
List rcpp_GridSamplePhi(NumericVector T, NumericVector u, NumericVector chi2, NumericVector tau, NumericVector e, NumericVector logl, NumericVector nu)
{
  //
  int i, j;
  int nProps = T.size();
  int n      = e.size();
  double phi;
  double v2;
  double Lprop_min = 1e100;
  
  NumericVector Lprop(nProps);
  
  // Compute negative log-posteriors for each proposal
  for (i = 0; i<nProps; i++)
  {
    Lprop[i] = 0;
    // Get proposed parameter values
    phi = T(i);
    
    for (j = 0; j<n; j++)
    {
      v2 = chi2[0]*(phi*phi + (1-phi)*(1-phi)*exp(logl[j]*2.0*tau[0]));
      Lprop[i] = Lprop[i] + (nu[0]+1.0)*0.5*log(1 + e[j]*e[j]/v2/nu[0]) + 0.5*log(v2);
      //Lprop[i] = Lprop[i] + 0.5*e[j]*e[j]/v2 + log(v2)*0.5;
    }
    if (Lprop[i] < Lprop_min)
    {
      Lprop_min = Lprop[i];
    }
  }
  
  // Sample
  i = grid_sample(u[0], Lprop, Lprop_min);
  List rv = List::create(Named("theta") = T(i), 
                         Named("theta.ix") = i+1
  );
  
  return rv;
}



// [[Rcpp::export]]
List rcpp_GridSampleRho(NumericVector R, NumericVector u, NumericVector ytilde, NumericVector v2, NumericVector logl, NumericVector w1, NumericVector nu, NumericVector rhoscale)
{
  // rho.prop, runif(1), v2, y, l, b, log.l, w)
  int i, j;
  int nProps = R.size();
  int n      = ytilde.size();
  double Lprop_min = 1e100;
  double e;

  NumericVector Lprop(nProps);
  
  // Compute negative log-posteriors for each proposal
  for (i = 0; i<nProps; i++)
  {
    Lprop[i] = log(R[i]*R[i] + 1.0);
    // Get proposed parameter values
    for (j = 0; j<n; j++)
    {
      //v2 = chi2[0]*(phi*phi + (1-phi)*(1-phi)*exp(logl[j]*2.0*tau))*omega2[j];
      //e = y[j] - l[j] - w[1]*exp(R(i)*logl[j]) - w[2]*b[j];
      e = ytilde[j] - w1[0]*rhoscale[i]*exp(R(i)*logl[j]);
      //Lprop[i] = Lprop[i] + 0.5*e*e/v2[j];
      //Lprop[i] = Lprop[i] + (nu[0]+1)/2*log(1 + e*e/chi2[0]/nu[0]);
      Lprop[i] = Lprop[i] + (nu[0]+1)/2*log(1 + e*e/v2[j]/nu[0]);
    }
    if (Lprop[i] < Lprop_min)
    {
      Lprop_min = Lprop[i];
    }
  }
  
  // Sample
  i = grid_sample(u[0], Lprop, Lprop_min);
  List rv = List::create(Named("theta") = R(i), 
                         Named("theta.ix") = i+1
  );
  
  return rv;
}

// [[Rcpp::export]]
List rcpp_GridSampleRhoGaussianMix(NumericVector R, NumericVector u, NumericVector ytilde, NumericVector v2, NumericVector logl, NumericVector w1)
{
  // rho.prop, runif(1), v2, y, l, b, log.l, w)
  int i, j;
  int nProps = R.size();
  int n      = ytilde.size();
  double Lprop_min = 1e100;
  double e;
  
  NumericVector Lprop(nProps);
  
  // Compute negative log-posteriors for each proposal
  for (i = 0; i<nProps; i++)
  {
    Lprop[i] = log(R[i]*R[i] + 1.0);
    // Get proposed parameter values
    for (j = 0; j<n; j++)
    {
      e = ytilde[j] - w1[0]*exp(R(i)*logl[j]);
      Lprop[i] = Lprop[i] + 0.5*e*e/v2[j];
    }
    if (Lprop[i] < Lprop_min)
    {
      Lprop_min = Lprop[i];
    }
  }
  
  // Sample
  i = grid_sample(u[0], Lprop, Lprop_min);
  List rv = List::create(Named("theta") = R(i), 
                         Named("theta.ix") = i+1
  );
  
  return rv;
}
  
// [[Rcpp::export]]
List rcpp_expsmooth(NumericVector y, NumericVector alphaV, NumericVector betaV, NumericVector l1, NumericVector b1)
{
  int n = y.size();

  double alpha = alphaV[0];
  double beta  = betaV[0];
  
  // Allocate variables
  NumericVector l(n);
  NumericVector b(n);
  NumericVector dl_dalpha(n);
  NumericVector db_dalpha(n);
  NumericVector db_dbeta(n);
  
  NumericVector dl_dl1(n);
  NumericVector db_dl1(n);
  NumericVector db_db1(n);
  
  // Initialise
  dl_dl1[0] = 1;
  db_dl1[0] = 0;
  db_db1[0] = 1;
  
  l1        = y[0];
  l[0]      = l1[0];
  b[0]      = b1[0];
  
  // Exponential smoothing step
  int i;
  for (i = 1; i<n; i++)
  {
    // Update the smoothed curves
    l[i] = alpha*y[i] + (1-alpha)*l[i-1];
    b[i] = beta*(l[i] - l[i-1]) + (1-beta)*b[i-1];
    
    // Update gradients
    dl_dalpha[i] = y[i] + (1-alpha)*dl_dalpha[i-1] - l[i-1];
    db_dalpha[i] = beta*(dl_dalpha[i] - dl_dalpha[i-1]) + (1-beta)*db_dalpha[i-1];
    db_dbeta[i]  = (l[i] - l[i-1]) + (1-beta)*db_dbeta[i-1] - b[i-1];
    
    dl_dl1[i]    = (1-alpha)*dl_dl1[i-1];
    db_dl1[i]    = beta*(dl_dl1[i] - dl_dl1[i-1]) + (1-beta)*db_dl1[i-1];
    
    db_db1[i]    = (1-beta)*db_db1[i-1];
  }
  
  // Return results as a list
  List rv = List::create(Named("l") = l, 
                         Named("b") = b, 
                         Named("dl_dalpha") = dl_dalpha,
                         Named("db_dalpha") = db_dalpha,
                         Named("db_dbeta") = db_dbeta,
                         Named("dl_dl1") = dl_dl1,
                         Named("db_dl1") = db_dl1,
                         Named("db_db1") = db_db1
                        );
  
  return rv;
}

// [[Rcpp::init]]
void rstan_additional_init(DllInfo *dll){
  R_useDynamicSymbols(dll, TRUE); // necessary for .onLoad() to work
}
