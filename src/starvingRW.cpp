#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List starvingRW(int L, int s_max, int s_crit, int gain, int t_max, double pr_grow, double pr_rep, double pr_mort, IntegerVector srw, IntegerVector rwloc, IntegerVector r) {
  
  IntegerVector pop_c(t_max);
  IntegerVector pop_r(t_max);
  
  //Begin time loop
  for (int t=0; t<t_max; t++) {
    
    //Across each individual in the system...
    bool ind_check = TRUE;
    int num = srw.size();
    pop_c[t] = num;
    pop_r[t] = sum(r);
    int i = -1;
    
    while (ind_check == TRUE) {
      
      //Individual index
      i = i + 1;
      
      //Move to a new location
      IntegerVector nn(4);
      int x = rwloc(i);
      int check = 0;
      //Bottom Row
      if (x <= L+2) {int new_x = x + L*(L+2); check = 1;}
      //Top Row
      if (x >= ((L+2)^2 - (L+1))) {int new_x = x - L*(L+2); check = 1;}
      //Left Row
      if (x%(L+2) == 1) {int new_x = x + L; check = 1;}
      //Right Row
      if (x%(L+2) == 0) {int new_x = x - L; check = 1;}
      //Bottom Left
      if ((x <= L+2) && (x%(L+2) == 1)) {int new_x = x + L*(L+2) + L; check = 1;}
      //Bottom Right
      if ((x <= L+2) && (x%(L+2) == 0)) {int new_x = x + L*(L+2) - L; check = 1;}
      //Top Left
      if ((x >= ((L+2)^2 - (L+1))) && (x%(L+2) == 1)) {int new_x = x - L*(L+2) + L; check = 1;}
      //Top Right
      if ((x >= ((L+2)^2 - (L+1))) && (x%(L+2) == 0)) {int new_x = x - L*(L+2) - L; check = 1;}
      //Middle
      if(check == 0) {int new_x = x}
      IntegerVector nn(4); nn(0) = new_x+1; nn(1) = new_x-1; nn(2) = new_x+(L+2); nn(3) = new_x-(L+2);
      //Sample new_x
      int num_draw = 4;
      NumericVector rdraw_v = runif(1);
      double rdraw = as<double>(rdraw_v);
      int draw = (int) floor(rdraw*(num_draw));
      int new_loc = nn(draw);
      
      //Consume resource if it is there
      int new_s = min(srw(i) - 1 + r(new_loc)*gain,s_max);
      //deplete the resource
      r(new_loc) = 0;
      
      //Impliment mortality and update
      //Random mortality
      if (new_s <= 0) {
        //Draw random value
        NumericVector rand = runif(1);
        double rdraw = as<double>(r_draw);
        if (rdraw < pr_mort) {
          //Spaghetti method of removal
          srw(i) = srw(num);
          IntegerVector new_srw(num-1);
          for (int k=0;k<(num-1);k++) {
            new_srw(k) = srw(k);
          }
          IntegerVector srw = new_srw;
          rwloc(i) = rwloc(num);
          IntegerVector new_rwloc(num-1);
          for (int k=0;k<(num-1);k++) {
            new_rwloc(k) = rwloc(k);
          }
          IntegerVector rwloc = new_rwloc;
          i = i - 1; //reanalyzes the 'new member' of the slot
        } else {
          //Individual survives
          //Record the new location and new state
          rwloc(i) = new_loc;
          srw(i) = new_s;
        }
      } else {
        //Individual survives:
        //Record the new location and new state
        rwloc(i) = new_loc;
        srw(i) = new_s;
      }
      
      //Recalculate num (number of individuals)... will be shorter if there is mortality
      int num = srw.size();
      if (i == num) {
        bool ind_check = FALSE //The while loop is stopped if you get to the end of the list
      }
      
    } //End While loop
    
    
    
    
    
    
  } //End t loop
  
  
  
}
