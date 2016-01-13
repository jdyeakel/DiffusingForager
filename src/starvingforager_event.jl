
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
DF = 0.2;
DH = 0.2;

#Starting number of foragers
initsize = 10;
#Initial state values
ind_vec = rand(collect(0:2),initsize);
#Initial location values
loc_vec = rand(collect(1:size),initsize);

#NEED TO ENSURE THAT RESOURCES ARE PLACED 1 PER SITE
#HAVE A noresourcesites vector that accounts for sites WITHOUT resources
#Updated each timestep


t = 0;

#Initial count of how many resouces, starvers, and full in this timestep??
#Count the number of individual R + S + F
tot = length(ind_vec);
NR = length(find(x->x==0,ind_vec));
NH = length(find(x->x==1,ind_vec));
NF = initsize - (R + S);


tic = 1;
while t < (t_term-1)

  N = NF + NH + NR;
  F = NF/N;
  H = NH/N;
  R = NR/N;


  #Calculate Rate
  Rate = F*(lambda + sigma*(1-R) + DF) + H*(rho*R + mu + DH) + R*(alpha*(1-R) + (F+H));
  dt = 1/(Rate*N);


  #Construct probability lines, which are a function of R, S, F

  #Grow <-----> Starve <-----> Diffuse
  F_pr_line = Array(Float64,2);
  F_pr_line[1] = lambda/(lambda+sigma*(K-R) + DF);
  F_Pr_line[2] = F_pr_line[1] + (sigma*(K-R))/(lambda+sigma*(K-R) + DF);

  #Recover <-----> Mortality <-----> Diffuse
  H_pr_line = Array(Float64,2);
  H_pr_line[1] = (rho*R)/(rho*R + mu + DH);
  H_pr_line[2] = mu/(rho*R + mu + DH);

  #Grow <-----> Consumed
  R_pr_line = (alpha*(K-R))/((alpha*(K-R)) + (F + H));



  #Randomly select an individual (R,S,F) with probability 1/N
  #ind thus represents the POSITION of the individual
  #Update total number of individuals

  #Update the total
  tot = length(ind_vec);
  id = rand(collect(1:tot))

  #If randomly selected individual is Full
  if ind_vec(id) == 2
    state = 2;
    location = loc_vec[id];

    draw_event = rand();

    #GROW
    if draw_event < F_pr_line[1]
      push!(ind_vec,state);
      push!(loc_vec,location);
      F = F+1;
    end

    #STARVE
    if draw_event >= F_pr_line[1] && draw_event < F_pr_line[2]
      ind_vec[id] = 1;

      S = S+1;
      F = F-1;
    end

    #MOVE
    if draw_event > F_pr_line[2]
      #Choose random new location
      loc_vec[ind] = rand(collect(1:size))
    end


  end



  #If randomly selected individual is Hungry
  if ind_vec(id) == 1
    state = 1;
    location = loc_vec[id];

    #Draw a random event
    #Recover, die, or move??
    draw_event = rand();

    #Recover
    if draw_event < H_pr_line[1]
      loc_vec[id] = 2;
      S = S-1;
      F = F+1;
    end

    #Die
    if draw_event >= H_pr_line[1] && draw_event < H_pr_line[2]
      deleteat!(ind_vec,id);
      deleteat!(loc_vec,id);
      S = S-1;
    end

    #MOVE
    if draw_event > H_pr_line[2]
      loc_vec[ind] = rand(collect(1:size))
    end

  end


  #If randomly selected individual is a resource
  if ind_vec(id) == 0
    state = 0;

    #Draw a random event
    #Grow, become consumed or move?
    draw_event = rand();

    #GROW
    if draw_event < R_Pr_line

      #Append a new resource to the END of the vector
      push!(ind_vec,state);

      #Choose random no-resource site
      newresourcepos = rand(collect(1:length(noresourcesites)));
      location = noresourcesites[newresourcepos];

      #Append the resource's location to the END of the vector
      push!(loc_vec,location);
      #Update Tally
      R = R + 1;

      #Update the no resource site by deleting the position that has been filled
      deleteat!(noresourcesites,newresourcepos);
      
    end

    #BECOME CONSUMED!
    if draw_event >= R_pr_line && draw_event < 1
      deleteat!(ind_vec,id);
      deleteat!(loc_vec,id);
      R = R-1;
    end


  end




  #NEED TO DEFINE dt

  #Advance time
  t = t + dt;

  #Update output



end #end while loop over t


#Posthoc analysis
