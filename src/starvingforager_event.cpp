#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List starvingforager_event(
  int L,
  int t_term,
  double sigma,
  double rho,
  double lambda,
  double mu
) {
  double t = 0;

  //ind_vec: the vector of individual states... 0 = resource, 1=starver, 2=full
  //pos_vec: the vector of individual locations

  while (t < t_term) {
    //Count the number of individual R + S + F
    int tot = ind_vec.size();

    //Randomly select an individual (R,S,F) with probability 1/N
    //ind thus represents the POSITION of the individual
    int id = runif(1,0,tot-1);

    //If ind is a resource...
    if (ind_vec(id) == 0) {

    }
    //If ind is a starver...
    if (ind_vec(id) == 1) {

    }
    //If ind is a consumer...
    if (ind_vec(id) == 2) {
      
    }





  }


}
