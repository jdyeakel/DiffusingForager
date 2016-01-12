
#Read in packages/function
#ipbc :: torus movement
include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/ipbc.jl")


#Initiate values
#Lattice dimension
L = 20;
dim = 2;
size = (L-2)^dim;
t_term = 500;

alpha = 0.5;
K = 1;
sigma = 0.25;
rho = 0.25;
lambda = 0.3;
mu = 0.1;

#Starting number of foragers
initsize = 10;
#Initial state values
ind_vec = rand(collect(0:2),initsize);
#Initial location values
loc_vec = rand(collect(1:size),initsize);


t = 0;

#Initial count of how many resouces, starvers, and full in this timestep??
#Count the number of individual R + S + F
tot = length(ind_vec);
R = length(find(x->x==0,ind_vec));
S = length(find(x->x==1,ind_vec));
F = initsize - (R + S);


tic = 1;
while t < (t_term-1)

  #Construct probability lines, which are a function of R, S, F

  #Grow <-----> Consumed
  R_pr_line = (alpha*(K-R))/((alpha*(K-R)) + (F + S));

  #Recover <-----> Mortality
  S_pr_line = (rho*R)/(rho*R + mu);

  #Grow <-----> Starve
  F_pr_line = lambda/(lambda+sigma*(K-R));

  #Randomly select an individual (R,S,F) with probability 1/N
  #ind thus represents the POSITION of the individual
  #Update total number of individuals

  #Update the total
  tot = length(ind_vec);
  id = rand(collect(1:tot))

  #If randomly selected individual is a resource
  if ind_vec(id) == 0
    #Append a new resource to the END of the vector
    push!(ind_vec,state);
    #Append the resource's location to the END of the vector
    push!(loc_vec,location);
    #Update Tally
    R = R + 1;

  end

  #If randomly selected individual is a starver
  if ind_vec(id) == 1


  end

  #If randomly selected individual is a full
  if ind_vec(id) == 2


  end
