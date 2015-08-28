#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List starvingRW_pr(
  int L,
  int s_max,
  int s_crit,
  int gain,
  int t_max,
  double alpha,
  double sigma,
  double rho,
  double lambda,
  double mu,
  IntegerVector srw,
  IntegerVector rwloc,
  IntegerVector r,
  double p) {

  //In this simulation, if p=1 it is in the continuous state
  //If p=0, we are in the fully spatial state
  //Intermediate values mean that with probability p, you function in
  //accordance to the mean field model


  IntegerVector pop_c(t_max);
  IntegerVector pop_r(t_max);
  IntegerVector pop_full(t_max);
  IntegerVector pop_starve(t_max);
  //List r_frame(t_max);
  int rsize = r.size();
  IntegerMatrix rm(t_max,rsize);
  IntegerMatrix locm(t_max,rsize);

  int L_size = pow((L+2),2);

  //Begin time loop
  for (int t=0; t<t_max; t++) {

    //Across each individual in the system...
    int ind_check = 1;
    int num = srw.size();
    pop_c(t) = num;
    pop_full(t) = sum(srw);
    pop_starve(t) = num - pop_full(t);
    pop_r(t) = sum(r);
    //Rcout << "got here! t... " << r(100) << std::endl;
    for (int j=0;j<rsize;j++) {
      rm(t,j) = r(j);
    }
    for (int j=0;j<num;j++) {
      int l = rwloc(j);
      locm(t,l) = 1;
    }
    //r_frame(t) = r;

    //Resource proportions at time t
    double total_r = (double) sum(r);
    double total_sites = (double) L_size;
    double pRr = total_r / total_sites; //Proportion of sites with resource
    double pRo = 1.0 - pRr; //Proportion of sites without resource

    //Define the probability of starvation and probability of recovery
    double pr_starve = sigma*pRo;
    double pr_recover = rho*pRr;
    double pr_grow = alpha*pRr;

    //Rcout << "Prob(resource growth) = " << pr_starve << std::endl;

    int i = 0;
    int x;
    while (ind_check == 1) {

      NumericVector pdraw_v = runif(1);
      double pdraw = as<double>(pdraw_v);

      //Rcout << "beginningwhile =  " << i << std::endl;
      //Rcout << "rwlocsize...  " << rwloc.size() << std::endl;
      //Rcout << "srwsize...  " << srw.size() << std::endl;
      //Rcout << "num...  " << num << std::endl;

      //With p probability, move to a random location
      //With (1-p) probability, move to a nn location

      int new_loc;
      if (pdraw < p) {
        NumericVector sitedraw_v = runif(1);
        double sitedraw = as<double>(sitedraw_v);
        int sitedraw_rv = (int) floor(sitedraw*(L_size));
        new_loc = sitedraw_rv;
      } else {
        //With (1-p) probability, move to a nn location

        //Rcout << "got here! 1... " << t << std::endl;
        //Move to a new location
        IntegerVector nn(4);
        //Rcout << "Mortality! 1 " << ind_check << std::endl;
        x = rwloc(i);
        //Convert to {1:L}
        x = x+1;
        int new_x;
        int check = 0;
        int mod = (x)%((L+2));
        //Bottom Row
        if (x <= (L+2)) {new_x = x + L*(L+2); check = 1;}
        //Top Row
        if (x >= (pow((L+2),2) - (L+1))) {new_x = x - L*(L+2); check = 1;}
        //Left Row
        if (mod == 1) {new_x = x + L; check = 1;}
        //Right Row
        if (mod == 0) {new_x = x - L; check = 1;}
        //Bottom Left
        if ((x <= (L+2)) && (mod == 1)) {new_x = x + L*(L+2) + L; check = 1;}
        //Bottom Right
        if ((x <= (L+2)) && (mod == 0)) {new_x = x + L*(L+2) - L; check = 1;}
        //Top Left
        if ((x >= (pow((L+2),2) - (L+1))) && (mod == 1)) {new_x = x - L*(L+2) + L; check = 1;}
        //Top Right
        if ((x >= (pow((L+2),2) - (L+1))) && (mod == 0)) {new_x = x - L*(L+2) - L; check = 1;}
        //Middle
        if(check == 0) {new_x = x;}
        //Convert back to {0:L-1}
        new_x = new_x - 1;
        nn(0) = new_x+1; nn(1) = new_x-1; nn(2) = new_x+(L+2); nn(3) = new_x-(L+2);
        //Sample new_x
        int num_draw = 4;
        NumericVector rdraw_v = runif(1);
        double rdraw = as<double>(rdraw_v);
        int draw = (int) floor(rdraw*(num_draw));
        new_loc = nn(draw);
      }
      //Rcout << "I got here! t= " << t << std::endl;
      //Consume resource if it is there
      //Make this probabilistic... with probability sigma/rho (!!!!!!!!!!)
      int new_s;
      if (r(new_loc) == 1) {
        //If you find a resource, and are Full
        if (srw(i) == 1) {
          new_s = 1;
        }
        //If you find a resource, and are Starving, recover with rate rho
        if (srw(i) <= 0) {
          NumericVector rdraw_v = runif(1);
          double rdraw = as<double>(rdraw_v);
          if (rdraw < pr_recover) {
            new_s = 1;
          } else {
            new_s = 0;
          }
        }
      } else {
        //If you DONT find a resource, and are Full, starve with rate sigma
        if (srw(i) == 1) {
          NumericVector rdraw_v = runif(1);
          double rdraw = as<double>(rdraw_v);
          if (rdraw < pr_starve) {
            new_s = 0;
          } else {
            new_s = 1;
          }
        }
        //If you DONT find a resource and are starving, no change
        if (srw(i) <= 0) {
            new_s = 0;
        }
      }
      //IntegerVector t_vec(2);
      //t_vec(0) = srw(i) - 1 + r(new_loc)*gain;
      //t_vec(1) = s_max;
      //IntegerVector::iterator it = std::min_element(t_vec.begin(), t_vec.end());
      //int new_s = *it;
      //deplete the resource
      r(new_loc) = 0;

      //Impliment mortality and update
      //Random mortality
      if (new_s <= 0) {
        //Draw random value
        NumericVector rand = runif(1);
        double rdraw = as<double>(rand);
        if (rdraw < mu) {
          //Mortality occurs...


          //Rcout << "num =  " << num << std::endl;
          //erase the ith element of the srw vector
          srw.erase(srw.begin() + i);
          //erase the ith element of the rwloc vector
          rwloc.erase(rwloc.begin() + i);

//          //Spaghetti method of removal
//          int srw_last = srw(num-1);
//          srw(i) = srw_last; //this is the last element of the vector... because vector starts with 0
//          IntegerVector new_srw(num-1);
//          for (int k=0;k<(num-1);k++) {
//            new_srw(k) = srw(k);
//          }
//
//          IntegerVector srw = new_srw;
//          int rwloc_last = rwloc(num-1);
//          rwloc(i) = rwloc_last; //this is the last element of the vector... because vector starts with 0
//          IntegerVector new_rwloc(num-1);
//          for (int k=0;k<(num-1);k++) {
//            new_rwloc(k) = rwloc(k);
//          }
//
//          IntegerVector rwloc = new_rwloc;
//          i = i - 1; //reanalyzes the 'new member' of the slot

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
      num = srw.size();

      //Individual index
      i = i + 1;

      //Check while loop
      //Make sure this is working how it is supposed to...
      if (i >= (num-1)) {
        ind_check = 0; //The while loop is stopped if you get to the end of the list
      }

    } //End While loop

    //4) RW reproduction
    for(int i=0;i<num;i++) {

      NumericVector pdraw_v = runif(1);
      double pdraw = as<double>(pdraw_v);

      int state = srw(i);
      if (state > s_crit) {
        NumericVector rand = runif(1);
        double rdraw = as<double>(rand);
        if (rdraw < lambda) {
          int srw_new = srw(i);
          int size_new = srw.size()+1;

          //With probability p seed offspring to a random site
          int rwloc_new;
          if (pdraw < p) {
            NumericVector offsitedraw_v = runif(1);
            double offsitedraw = as<double>(offsitedraw_v);
            int offsite_rv = (int) floor(offsitedraw*(L_size));
            rwloc_new = offsite_rv;
          } else {
            //With probability 1-p, seed offspring to the same site as parent
            rwloc_new = rwloc(i);
          }

          IntegerVector srw_newvec(size_new);
          IntegerVector rwloc_newvec(size_new);
          for (int k=0;k<size_new;k++) {
            if (k==(size_new - 1)) {
              srw_newvec(k) = srw_new;
              rwloc_newvec(k) = rwloc_new;
            } else {
              srw_newvec(k) = srw(k);
              rwloc_newvec(k) = rwloc(k);
            }
          }

          //redefine srw and rwloc
          srw = srw_newvec;
          rwloc = rwloc_newvec;
          //Rcout << "Reproduction! " << srw_newvec.size() << std::endl;
        }

      }

    }

    //Recalculate the number of individuals... will be greater if there is reproduction
    num = srw.size();


    //Resource growth
    for(int j=0; j<L_size; j++) {
      //Only determine growth if r(j) = 0...

      if (r(j) == 0) {
        NumericVector pdraw_v = runif(1);
        double pdraw = as<double>(pdraw_v);

        int r_nn;
        if (pdraw < p) {

          NumericVector randv = runif(1);
          double rdraw = as<double>(randv);
          //The probability of growth is scaled to the total num resources...
          //pr_grow_scaled = 1 - pow((1 - pr_grow),num_r);

          if (rdraw < pr_grow) {
            r(j) = 1;
          }
        } else {
          //Determine the number of nearest neighbor resources
          IntegerVector nn(4);
          int x = j;
          //Convert to {1:L}
          x = x+1;
          int new_x;
          int check = 0;
          int mod = (x)%((L+2));
          //Bottom Row
          if (x <= (L+2)) {new_x = x + L*(L+2); check = 1;}
          //Top Row
          if (x >= (pow((L+2),2) - (L+1))) {new_x = x - L*(L+2); check = 1;}
          //Left Row
          if (mod == 1) {new_x = x + L; check = 1;}
          //Right Row
          if (mod == 0) {new_x = x - L; check = 1;}
          //Bottom Left
          if ((x <= (L+2)) && (mod == 1)) {new_x = x + L*(L+2) + L; check = 1;}
          //Bottom Right
          if ((x <= (L+2)) && (mod == 0)) {new_x = x + L*(L+2) - L; check = 1;}
          //Top Left
          if ((x >= (pow((L+2),2) - (L+1))) && (mod == 1)) {new_x = x - L*(L+2) + L; check = 1;}
          //Top Right
          if ((x >= (pow((L+2),2) - (L+1))) && (mod == 0)) {new_x = x - L*(L+2) - L; check = 1;}
          //Middle
          if(check == 0) {new_x = x;}
          //Convert back to {0:L-1}
          new_x = new_x - 1;
          nn(0) = new_x+1; nn(1) = new_x-1; nn(2) = new_x+(L+2); nn(3) = new_x-(L+2);
          //Total number of nearest neighbors to r
          r_nn = r(nn(0)) + r(nn(1)) + r(nn(2)) + r(nn(3));

          //The probability of growth is scaled to the total num resources...
          //pr_grow_scaled = 1 - pow((1 - pr_grow),r_nn);

          if (r_nn >= 1) {
            NumericVector rand = runif(1);
            double rdraw = as<double>(rand);
            if (rdraw < pr_grow) {
              r(j) = 1;
              //Rcout << "Mortality! 2 " << rdraw << std::endl;
            }
          }
        }
      } //end if
    } //End j

    //Rcout << "Mortality! 2 " << t << std::endl;

  } //End t loop

  List cout(7);
  cout(0) = pop_r;
  cout(1) = pop_c;
  cout(2) = pop_starve;
  cout(3) = pop_full;
  cout(4) = rm;
  cout(5) = locm;
  cout(6) = srw;

  return cout; //check

}

//ToDo: need to export Starving + Full consumers separately
//Build in mean field randomizations (separate file??)
