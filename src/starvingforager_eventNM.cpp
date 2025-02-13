#include <Rcpp.h>
using namespace Rcpp;

// NOTES: make proportions relative to lattice size
// Probably don't need count of RSF... just proportions

// [[Rcpp::export]]
List starvingforager_eventNM(
int L,          //Lattice dim
int t_term,     //Terminal time
double alpha,   //Resource growth rate
double K,       //Resource carrying capacity
double sigma,   //Starvation rate
double rho,     //Recovery rate
double lambda,  //Growth rate
double mu,      //Mortality rate
IntegerVector ind_vec, //Initial vector of states
IntegerVector loc_vec //Initial vector of locations
) {
    //Dimension of the lattice
    int dim = 2;
    //Lattice size
    double size = pow(L-2,dim);
    //Initial time
    double t = 0;

    double max;
    double min;


    //Output Lists
    List ind_out(1);
    List loc_out(1);
    NumericVector t_out(1);
    //The initial state
    ind_out(0) = ind_vec;
    loc_out(0) = loc_vec;
    t_out(0) = 0;
    //ind_vec: the vector of individual states... 0 = resource, 1=starver, 2=full
    //pos_vec: the vector of individual locations

    //Initial count of how many resouces, starvers, and full in this timestep??
    //Count the number of individual R + S + F
    int tot = ind_vec.size();
    double R = 0.L;
    double S = 0.L;
    double F = 0.L;
    // double Rp;
    // double Sp;
    // double Fp;
    for (int i=0;i<tot;i++) {
        if (ind_vec(i) == 0) {
            R = R + 1.L;
        }
        if (ind_vec(i) == 1) {
            S = S + 1.L;
        }
        if (ind_vec(i) == 2) {
            F = F + 1.L;
        }
    }

    //R,S,F are thus densities over the landscape of size 'size'
    // R = R/size;
    // S = S/size;
    // F = F/size;

    double R_pr_line;
    double S_pr_line;
    double F_pr_line;

    //Iterate over time
    //The loop stops when t > t_term-1... and will record the last value
    int tic = 1;
    while (t < (t_term-1)) {

        //Construct probability lines, which are a function of R, S, F

        //Grow <-----> Consumed
        R_pr_line = (alpha*(K-R))/((alpha*(K-R)) + (F + S));
        //R_pr_line(1) = R_pr_line(0) + ((F + S)/((alpha*(K-R)) + (F + S) + Dr));
        //R_pr_line(2) = R_pr_line(1) + (Dr/((alpha*(K-R)) + (F + S) + Dr));

        //Recover <-----> Mortality
        S_pr_line = (rho*R)/(rho*R + mu);
        //S_pr_line(1) = S_pr_line(0) + (mu/(rho*R + mu + Ds));
        //S_pr_line(2) = S_pr_line(1) + (Ds/(rho*R + mu + Ds));

        //Grow <-----> Starve
        F_pr_line = lambda/(lambda+sigma*(K-R));
        //F_pr_line(1) = F_pr_line(0) + ((sigma*(K-R))/(lambda+sigma*(K-R)+Df));
        //F_pr_line(2) = F_pr_line(1) + (Df/(lambda+sigma*(K-R)+Df));


        //Initiate variables
        double dt;

        //Randomly select an individual (R,S,F) with probability 1/N
        //ind thus represents the POSITION of the individual
        //Update total number of individuals
        tot = ind_vec.size();
        max = (double)(tot - 1);
        min = 0.L;
        int id = min + (rand() % (int)(max - min + 1));

        int state;
        int location;

        double draw_event;
        //If ind is a resource...
        if (ind_vec(id) == 0) {
            state = 0;
            location = loc_vec(id);

            //Draw a random event
            //Grow, become consumed or move?
            draw_event = ((double) rand() / (RAND_MAX));

            //Grow
            if (draw_event < R_pr_line) {
                //Append a new resource to the END of the vector
                ind_vec.push_back(state);
                //Append the resource's location to the END of the vector
                loc_vec.push_back(location);
                //Update Tally
                R = R + 1; //(1.L/size);
            }
            //Become consumed!!!!
            if ((draw_event >= R_pr_line) && (draw_event < 1.L)) {
                //Remove the consumed resource from the state vector
                ind_vec.erase(id);
                //Remove the consumed resource form the location vector
                loc_vec.erase(id);
                //Update Tally
                R = R - 1; //(1.L/size);
            }
            dt = 1.L/((alpha*(K-R)) + (F + S));
        }
        //If ind is a starver...
        if (ind_vec(id) == 1) {
            state = 1;
            location = loc_vec(id);

            //Draw a random event
            //Recover, die, or move??
            draw_event = ((double) rand() / (RAND_MAX));

            //Recover
            if (draw_event < S_pr_line) {
                //Update the state from starver to full
                ind_vec(id) = 2;
                //Update Tally
                S = S - 1; //(1.L/size);
                F = F + 1; //(1.L/size);
            }
            //Die
            if ((draw_event >= S_pr_line) && (draw_event < 1.L)) {
                //Remove the consumed resource from the state vector
                ind_vec.erase(id);
                //Remove the consumed resource form the location vector
                loc_vec.erase(id);
                //Update Tally
                S = S - 1; //(1.L/size);
            }
            dt = 1.L/(rho*R + mu);
        }
        //If ind is Full...
        if (ind_vec(id) == 2) {
            state = 2;
            location = loc_vec(id);

            //Draw a random event
            //Grow, starve, or move?
            draw_event = ((double) rand() / (RAND_MAX));

            //Grow
            if (draw_event < F_pr_line) {
                //Append a new resource to the END of the vector
                ind_vec.push_back(state);
                //Append the resource's location to the END of the vector
                loc_vec.push_back(location);
                F = F + 1; //(1.L/size);
            }
            //Starve
            if ((draw_event >= F_pr_line) && (draw_event < 1.L)) {
                //Update the state from full to starver
                ind_vec(id) = 1;
                //Update Tally
                F = F - 1; //(1.L/size);
                S = S + 1; //(1.L/size);
            }
            dt = 1.L/(lambda+sigma*(K-R));
        }

        //Advance time
        t = t + dt;
        //Rcout << "t = " << dt << std::endl;
        //Update output
        ind_out.push_back(ind_vec);
        loc_out.push_back(loc_vec);
        t_out.push_back(t);
        tic = tic + 1;

    } //end while loop over t

    List cout(3);
    cout(0) = ind_out;
    cout(1) = loc_out;
    cout(2) = t_out;
    return(cout);

}
