
#Read in packages/function
#ipbc :: torus movement
using StatsBase
using Gadfly
#include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/ipbc.jl")


#Initiate values
#Lattice dimension
L = 10;
dim = 2;
size = (L-2)^dim;
t_term = 50;

alpha = 0.5;
K = 1;
sigma = 0.25;
rho = 0.25;
lambda = 0.3;
mu = 0.1;
DF = 0.2;
DH = 0.2;

#Starting number of foragers & resources
initsize = 10;
#Initial state values
ind_vec = [rand(collect(1:2),Int(2*round(initsize/3)));zeros(Int64,Int(round(initsize/3)))];
#Initial location values
#2/3 of initsize will be foragers; 1/3 will be resources
resourceloc_vec = sample(collect(1:size),Int(round(initsize/3)),replace=false);
loc_vec = [rand(collect(1:size),Int(2*round(initsize/3)));resourceloc_vec];

#Re-establish initsize to account for rounding errors
initsize = length(ind_vec);

#Arrays for output
ind_out = Array(Array{Int64},1);
loc_out = Array(Array{Int64},1);
time_out = Array(Float64,1);

ind_out[1] = ind_vec;
loc_out[1] = loc_vec;
time_out[1] = 0;

#NEED TO ENSURE THAT RESOURCES ARE PLACED 1 PER SITE
#HAVE A noresourcesites vector that accounts for sites WITHOUT resources
#Updated each timestep
noresourcesites = collect(1:size);
deleteat!(noresourcesites,sort(resourceloc_vec));

t = 0;

#Initial count of how many resouces, starvers, and full in this timestep??
#Count the number of individual R + S + F
#tot = length(ind_vec);
NF = length(find(x->x==2,ind_vec));
NH = length(find(x->x==1,ind_vec));
NR = initsize - (NF + NH);
N = NF + NH + NR;

#Initial densities
F = NF/N;
H = NH/N;
R = NR/N;

tic = 0;

F_pr_line = Array(Float64,2);
H_pr_line = Array(Float64,2);
R_pr_line = Array(Float64,1);

id = Array(Float64,1);

while t < (t_term-1)
  tic = tic + 1;


  #Construct probability lines, which are a function of R, S, F
  #Grow <-----> Starve <-----> Diffuse

  F_pr_line[1] = lambda/(lambda+sigma*(K-R) + DF);
  F_pr_line[2] = F_pr_line[1] + (sigma*(K-R))/(lambda+sigma*(K-R) + DF);

  #Recover <-----> Mortality <-----> Diffuse

  H_pr_line[1] = (rho*R)/(rho*R + mu + DH);
  H_pr_line[2] = mu/(rho*R + mu + DH);

  #Grow <-----> Consumed
  R_pr_line = (alpha*(K-R))/((alpha*(K-R)) + (F + H));



  #Randomly select an individual (R,S,F) with probability 1/N
  #ind thus represents the POSITION of the individual
  #Update total number of individuals

  #Update the total

  #Choose a random individual position
  id = rand(collect(1:N));
  #Also, we need to make an ind_vec_old that is NOT updated for the sequence of IF statements
  ind_vec_old = ind_vec;

  #If randomly selected individual is Full
  if ind_vec_old[id] == 2
    state = 2;

    draw_event = rand();

    #GROW
    if draw_event < F_pr_line[1]
      push!(ind_vec,2); #ensure the new individual is in the full state
      #Offspring appears at a random location
      location = rand(collect(1:size));
      push!(loc_vec,location);
      NF = NF+1;
    end

    #STARVE
    if draw_event >= F_pr_line[1] && draw_event < F_pr_line[2]
      ind_vec[id] = 1;

      NH = NH+1;
      NF = NF-1;
    end

    #MOVE
    if draw_event > F_pr_line[2]
      #Move to a random location
      loc_vec[id] = rand(collect(1:size));
    end


  end



  #If randomly selected individual is Hungry
  if ind_vec_old[id] == 1
    state = 1;

    #Draw a random event
    #Recover, die, or move??
    draw_event = rand();

    #Recover (become full)
    if draw_event < H_pr_line[1]
      loc_vec[id] = 2;
      NH = NH-1;
      NF = NF+1;
    end

    #Die
    if draw_event >= H_pr_line[1] && draw_event < H_pr_line[2]
      deleteat!(ind_vec,id);
      deleteat!(loc_vec,id);
      NH = NH-1;
    end

    #MOVE
    if draw_event > H_pr_line[2]
      #Move to a random location
      loc_vec[id] = rand(collect(1:size));
    end

  end


  #If randomly selected individual is a resource
  if ind_vec_old[id] == 0
    state = 0;

    #Draw a random event
    #Grow, become consumed or move?
    draw_event = rand();

    #GROW
    if draw_event < R_pr_line

      #Append a new resource to the END of the vector
      push!(ind_vec,state);

      #Choose random no-resource site
      newresourcepos = rand(collect(1:length(noresourcesites)));
      location = noresourcesites[newresourcepos];

      #Append the resource's location to the END of the vector
      push!(loc_vec,location);
      #Update Tally
      NR = NR + 1;

      #Update the no resource site by deleting the position that has been filled
      deleteat!(noresourcesites,newresourcepos);

    else

      #BECOME CONSUMED!
      deleteat!(ind_vec,id);
      deleteat!(loc_vec,id);
      NR = NR-1;

    end

  end

  #Recalculate the size of the foragers and resources
  N = NF + NH + NR;
  #Recalculate the densities of each
  F = NF/N;
  H = NH/N;
  R = NR/N;

  #Calculate Rate
  Rate = F*(lambda + sigma*(1-R) + DF) + H*(rho*R + mu + DH) + R*(alpha*(1-R) + (F+H));
  dt = 1/(Rate*N);

  #Advance time
  t = t + dt;

  #Update output
  push!(ind_out,ind_vec);
  push!(loc_out,loc_vec);
  push!(time_out,t);

  #Break loop if extinction occurs
  if length(ind_vec) < 1
    print("Extinction has occured at t=",round(t,2)," and loop ",tic)
    break
  end

end #end while loop over t


#Posthoc analysis
steps = tic;
pop_F = zeros(steps);
pop_H = zeros(steps);
pop_R = zeros(steps);
for i=1:steps
  pop_F[i] = length(find(x->x==2,ind_out[i]));
  pop_H[i] = length(find(x->x==1,ind_out[i]));
  pop_R[i] = length(find(x->x==0,ind_out[i]));
end
