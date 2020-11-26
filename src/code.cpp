#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// [[Rcpp::export]]
Rcpp::List compute_cov(int i, int l1, int l2, int K, arma::vec meanE, arma::mat muhats, arma::mat deltainv, arma::vec eta, arma::mat theta, arma::mat zeta, arma::mat Ipq) {
  arma::vec covs(9);

  arma::mat temp = (meanE[i-1] / eta[i-1]) * (trans(muhats.col(l1-1)) * deltainv * deltainv * muhats.col(l2-1));
  covs[0] = temp(0,0);

  if (i == l2) {
    temp = (meanE[i-1] / eta[i-1]) * (trans(muhats.col(l1-1)) * deltainv * deltainv * muhats.col(i-1));
    covs[1] = temp(0,0);
  } else {
    covs[1] = 0;
  }

  if (i == l1) {
    temp = (meanE[l1-1] / eta[l1-1]) * (trans(muhats.col(l1-1)) * deltainv * deltainv * muhats.col(l2-1));
    covs[2] = temp(0,0);
  } else {
    covs[2] = 0;
  }

  if (l1 == l2) {
    temp = (meanE[l1-1] / eta[l1-1]) * (trans(muhats.col(i-1)) * deltainv * deltainv * muhats.col(i-1));
    covs[3] = temp(0,0);
  } else {
    covs[3] = 0;
  }

  temp = (meanE[i-1]) * (trans(muhats.col(l1-1)) * deltainv * deltainv * muhats.col(i-1) * trans(muhats.col(l2-1)) * deltainv * muhats.col(i-1));
  covs[4] = temp(0,0);

  temp = (meanE[l1-1]) * (trans(muhats.col(i-1)) * deltainv * deltainv * muhats.col(i-1) * trans(muhats.col(l2-1)) * deltainv * muhats.col(l1-1));
  covs[5] = temp(0,0);

  temp = (meanE[i-1]) * (trans(muhats.col(i-1)) * deltainv * muhats.col(l1-1) * trans(muhats.col(i-1)) * deltainv * deltainv * muhats.col(l2-1));
  covs[6] = temp(0,0);

  temp = (meanE[l2-1]) * (trans(muhats.col(l2-1)) * deltainv * muhats.col(l1-1) * trans(muhats.col(i-1)) * deltainv * deltainv * muhats.col(i-1));
  covs[7] = temp(0,0);

  float cov9 = 0;
  float psiil11 = 0;
  float psiil21 = 0;
  float sigma2l1l12 = 0;
  float sigma2l2l22 = 0;
  float sigma2il13 = 0;
  float sigma2il14 = 0;
  float sigma2il15 = 0;
  float sigma2il23 = 0;
  float sigma2il24 = 0;
  float sigma2il25 = 0;
  for (int r = 1; r <= K; r++) {
    temp = eta[r-1]*(theta(i-1,r-1)*(1-theta(i-1,r-1))+theta(l1-1,r-1)*(1-theta(l1-1,r-1)))*(trans(muhats.col(i-1))*deltainv*Ipq*deltainv*muhats.col(l1-1));
    psiil11 = psiil11 + temp(0,0);
    temp = eta[r-1]*(theta(i-1,r-1)*(1-theta(i-1,r-1))+theta(l2-1,r-1)*(1-theta(l2-1,r-1)))*(trans(muhats.col(i-1))*deltainv*Ipq*deltainv*muhats.col(l2-1));
    psiil21 = psiil21 + temp(0,0);
    sigma2l1l12 = sigma2l1l12 + eta[r-1]*theta(l1-1,r-1)*(1-theta(l1-1,r-1))*zeta(l1-1,r-1)*zeta(l1-1,r-1)*(1/eta[l1-1]-2*zeta(l1-1,l1-1));
    sigma2l2l22 = sigma2l2l22 + eta[r-1]*theta(l2-1,r-1)*(1-theta(l2-1,r-1))*zeta(l2-1,r-1)*zeta(l2-1,r-1)*(1/eta[l2-1]-2*zeta(l2-1,l2-1));
    sigma2il13 = sigma2il13 + eta[r-1]*theta(i-1,r-1)*(1-theta(i-1,r-1))*zeta(l1-1,r-1)*zeta(l1-1,r-1)*(1/eta[i-1]-2*zeta(i-1,i-1));
    sigma2il14 = sigma2il14 + eta[r-1]*theta(l1-1,r-1)*(1-theta(l1-1,r-1))*zeta(i-1,r-1)*zeta(i-1,r-1)*(1/eta[l1-1]-2*zeta(l1-1,l1-1));
    sigma2il15 = sigma2il15 + eta[r-1]*(theta(i-1,r-1)*(1-theta(i-1,r-1))+theta(l1-1,r-1)*(1-theta(l1-1,r-1)))*zeta(i-1,r-1)*zeta(r-1,l1-1)*zeta(i-1,l1-1);
    sigma2il23 = sigma2il23 + eta[r-1]*theta(i-1,r-1)*(1-theta(i-1,r-1))*zeta(l2-1,r-1)*zeta(l2-1,r-1)*(1/eta[i-1]-2*zeta(i-1,i-1));
    sigma2il24 = sigma2il24 + eta[r-1]*theta(l2-1,r-1)*(1-theta(l2-1,r-1))*zeta(i-1,r-1)*zeta(i-1,r-1)*(1/eta[l2-1]-2*zeta(l2-1,l2-1));
    sigma2il25 = sigma2il25 + eta[r-1]*(theta(i-1,r-1)*(1-theta(i-1,r-1))+theta(l2-1,r-1)*(1-theta(l2-1,r-1)))*zeta(i-1,r-1)*zeta(r-1,l2-1)*zeta(i-1,l2-1);
    temp = (eta[r-1]*meanE[r-1])*(trans(muhats.col(r-1))*deltainv*muhats.col(l1-1)*trans(muhats.col(i-1))*deltainv*deltainv*muhats.col(i-1)*trans(muhats.col(l2-1))*deltainv*muhats.col(r-1));
    cov9 = cov9 + temp(0,0);
  }

  covs[8] = cov9;

  float psiil12 = 0;
  float psiil22 = 0;
  float sigma2l1l13 = 0;
  float sigma2l2l23 = 0;
  float sigma2il16 = 0;
  float sigma2il26 = 0;
  for (int r = 1; r <= K; r++) {
    for (int s = 1; s <= K; s++) {
      temp = eta[r-1] * eta[s-1] * theta(s-1,r-1) * (1-theta(s-1,r-1)) * (trans(muhats.col(s-1))*deltainv*Ipq*deltainv*(muhats.col(l1-1)*trans(muhats.col(i-1))+muhats.col(i-1)*trans(muhats.col(l1-1)))*deltainv*muhats.col(s-1));
      psiil12 = psiil12 + temp(0,0);
      temp = eta[r-1] * eta[s-1] * theta(s-1,r-1) * (1-theta(s-1,r-1)) * (trans(muhats.col(s-1))*deltainv*Ipq*deltainv*(muhats.col(l2-1)*trans(muhats.col(i-1))+muhats.col(i-1)*trans(muhats.col(l2-1)))*deltainv*muhats.col(s-1));
      psiil22 = psiil22 + temp(0,0);
      sigma2l1l13 = sigma2l1l13 + eta[r-1]*eta[s-1]*theta(r-1,s-1)*(1-theta(r-1,s-1))*zeta(l1-1,r-1)*zeta(l1-1,r-1)*zeta(l1-1,s-1)*zeta(l1-1,s-1);
      sigma2l2l23 = sigma2l2l23 + eta[r-1]*eta[s-1]*theta(r-1,s-1)*(1-theta(r-1,s-1))*zeta(l2-1,r-1)*zeta(l2-1,r-1)*zeta(l2-1,s-1)*zeta(l2-1,s-1);
      sigma2il16 = sigma2il16 + eta[r-1]*eta[s-1]*theta(r-1,s-1)*(1-theta(r-1,s-1))*(zeta(i-1,r-1)*zeta(l1-1,s-1)+zeta(l1-1,r-1)*zeta(i-1,s-1))*(zeta(i-1,r-1)*zeta(l1-1,s-1)+zeta(l1-1,r-1)*zeta(i-1,s-1));
      sigma2il26 = sigma2il26 + eta[r-1]*eta[s-1]*theta(r-1,s-1)*(1-theta(r-1,s-1))*(zeta(i-1,r-1)*zeta(l2-1,s-1)+zeta(l2-1,r-1)*zeta(i-1,s-1))*(zeta(i-1,r-1)*zeta(l2-1,s-1)+zeta(l2-1,r-1)*zeta(i-1,s-1));
    }
  }

  float psiil1 = psiil11 - psiil12;
  float psiil2 = psiil21 - psiil22;
  float sigma2il1 = (theta(i-1,i-1)*(1-theta(i-1,i-1))+theta(l1-1,l1-1)*(1-theta(l1-1,l1-1)))*zeta(i-1,l1-1)*zeta(i-1,l1-1) + 2*theta(i-1,l1-1)*(1-theta(i-1,l1-1))*zeta(i-1,i-1)*zeta(l1-1,l1-1) + sigma2il13 + sigma2il14 - 2*sigma2il15 + sigma2il16/2;
  float sigma2il2 = (theta(i-1,i-1)*(1-theta(i-1,i-1))+theta(l2-1,l2-1)*(1-theta(l2-1,l2-1)))*zeta(i-1,l2-1)*zeta(i-1,l2-1) + 2*theta(i-1,l2-1)*(1-theta(i-1,l2-1))*zeta(i-1,i-1)*zeta(l2-1,l2-1) + sigma2il23 + sigma2il24 - 2*sigma2il25 + sigma2il26/2;
  float sigma2l1l1 = 4*theta(l1-1,l1-1)*(1-theta(l1-1,l1-1))*zeta(l1-1,l1-1)*zeta(l1-1,l1-1) + 4*sigma2l1l12 + 2*sigma2l1l13;
  float sigma2l2l2 = 4*theta(l2-1,l2-1)*(1-theta(l2-1,l2-1))*zeta(l2-1,l2-1)*zeta(l2-1,l2-1) + 4*sigma2l2l22 + 2*sigma2l2l23;

  arma::vec temps = {psiil1, psiil2, sigma2il1, sigma2il2, sigma2l1l1, sigma2l2l2};

  return Rcpp::List::create(covs, temps);
}







