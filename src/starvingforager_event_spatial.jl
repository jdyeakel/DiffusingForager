function starvingforager_event_spatial(L,dim,initsize,t_term,alpha,K,sigma,rho,lambda,mu,DF,DH)
  #Read in packages/function
  #ipbc :: torus movement
  include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/ipbc.jl")

  #Initiate values
  #Lattice dimension

  size = (L-2)^dim;


  #Initial state values
  ind_vec = [rand(collect(1:2),Int(2*round(initsize/3)));zeros(Int64,Int(round(initsize/3)))];
  #OR Pure Resource Start (for testing)
  #ind_vec = zeros(Int64,initsize);


  #Initial location values
  #2/3 of initsize will be foragers; 1/3 will be resources
  #replace = false because want only one resource per site
  resourceloc_vec = sample(collect(1:size),Int(round(initsize/3)),replace=false);
  loc_vec = [rand(collect(1:size),Int(2*round(initsize/3)));resourceloc_vec];

  #Re-establish initsize to account for rounding errors
  initsize = length(ind_vec);

  #Arrays for output (start out empty)
  # ind_out = (Array{Int64,1})[];
  # loc_out = (Array{Int64,1})[];
  time_out = Array(Float64,0);
  N_out = Array(Int64,0);
  #prop_out = (Array{Float64,1})[];

  #Copy so that ind_vec is independent of ind_out
  # push!(ind_out, copy(ind_vec));
  # push!(loc_out, copy(loc_vec));
  push!(time_out, 0);

  #NEED TO ENSURE THAT RESOURCES ARE PLACED 1 PER SITE
  #HAVE A noresourcesites vector that accounts for sites WITHOUT resources
  #Updated each timestep
  # noresourcesites = collect(1:size);
  # deleteat!(noresourcesites,sort(resourceloc_vec));

  t = 0;
  next_time_int = t + 1;

  #Initial count of how many resouces, starvers, and full in this timestep??
  #Count the number of individual R + S + F
  #tot = length(ind_vec);
  NF = length(find(x->x==2,ind_vec));
  NH = length(find(x->x==1,ind_vec));
  NR = initsize - (NF + NH);
  N = NF + NH + NR;

  #Initial densities
  F = NF/size;
  H = NH/size;
  R = NR/size;
  prop_out = Array{Float64}(3,1);
  prop_out[:,1] = [F,H,R];
  push!(N_out,N);
  #push!(prop_out,copy(prop));

  tic = 0;

  F_pr_line = Array(Float64,2);
  H_pr_line = Array(Float64,2);
  R_pr_line = Array(Float64,1);

  id = Array(Float64,1);

  while t < (t_term-1)

    #Note: we need to choose the individual first, because it is now the local density that determines the probabilities and the Rate (???)

    #Choose a random individual position
    id = rand(collect(1:N));

    #Also, we need to make an ind_vec_old that is NOT updated for the sequence of IF statements
    ind_vec_old = copy(ind_vec);

    #Nearest neighbor information
    #What is the current location?
    current_loc = loc_vec[id];
    #nearest neighbor sites
    nn_sites = ipbc(current_loc,L);
    l_nn = length(nn_sites);
    nn_states = zeros(Int64,4);
    for i=1:l_nn
      #what is the position of the individual at site nn_sites[i]?
      pos = find(x->x==nn_sites[i],loc_vec);
      #If the individual has a state...
      if length(pos) > 0
        nn_states[i] = ind_vec[pos[1]];
      #If the site is completely empty (no R, no H, no F), we will say it is '-1'
      else
        nn_states[i] = -1;
      end
    end

    #Local resource density ~ assumes 4 nearest neighbors
    R_locations = find(x->x==0,nn_states);
    R_local = length(R_locations)/l_nn;
    #What sites DON"T have resources?
    noresourcesites = nn_sites;
    deleteat!(noresourcesites,R_locations);

    #Calculate Rate
    Rate = F*(lambda + sigma*(1-R_local) + DF) + H*(rho*R_local + mu + DH) + R*(alpha*(1-R_local) + (F+H));
    dt = 1/(Rate*N);

    if t > next_time_int
      println("time= ",round(t,2))
      next_time_int = round(t+1,0)
    end


    if Rate > 0
      tic = tic + 1;

      #Construct probability lines, which are a function of R, S, F
      #Grow <-----> Starve <-----> Diffuse

      F_pr_line[1] = lambda/Rate;
      F_pr_line[2] = F_pr_line[1] + (sigma*(K-R_local))/Rate;

      #Recover <-----> Mortality <-----> Diffuse

      H_pr_line[1] = (rho*R_local)/Rate;
      H_pr_line[2] = H_pr_line[1] + mu/Rate;

      #Grow <-----> Consumed
      R_pr_line = (alpha*(K-R_local))/Rate;



      #Randomly select an individual (R,S,F) with probability 1/N
      #ind thus represents the POSITION of the individual
      #Update total number of individuals

      #Update the total




      #If randomly selected individual is Full
      if ind_vec_old[id] == 2
        state = 2;

        draw_event = rand();

        #GROW
        if draw_event < F_pr_line[1]
          push!(ind_vec,2); #ensure the new individual is in the full state
          #Offspring appears at SAME SITE
          location = loc_vec[id];
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
        if draw_event >= F_pr_line[2]
          #Move to a NEIGHBORING LOCATION
          loc_vec[id] = rand(nn_sites);
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
          ind_vec[id] = 2;
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
        if draw_event >= H_pr_line[2]
          #Move to a NEIGHBORING LOCATION
          loc_vec[id] = rand(nn_sites);
        end

      end


      #If randomly selected individual is a resource
      if ind_vec_old[id] == 0
        state = 0;

        #Draw a random event
        #Grow, become consumed or move?
        draw_event = rand();

        #GROW
        if draw_event < R_pr_line[1]

          #Append a new resource to the END of the vector
          push!(ind_vec,state);

          #Choose noresourcesite from the NEAREST NEIGHBORS
          location = rand(noresourcesites);

          #Append the resource's location to the END of the vector
          push!(loc_vec,location);
          #Update Tally
          NR = NR + 1;


        end

        #BECOME CONSUMED!
        if draw_event >= R_pr_line[1]

          #Delete individual from the list
          deleteat!(ind_vec,id);
          #Delete individual location from the list
          deleteat!(loc_vec,id);
          NR = NR-1;

        end

      end

      #Recalculate the size of the foragers and resources
      N = NF + NH + NR;
      #Recalculate the GLOBAL densities of each
      F = NF/size;
      H = NH/size;
      R = NR/size;

      prop = [F,H,R];

      #push!(prop_out,copy(prop));
      prop_out = hcat(prop_out,prop);

      #Advance time
      t = t + dt;



      #Make sure only values are pushed to the output
      ind_vec_new = copy(ind_vec);
      loc_vec_new = copy(loc_vec);

      #Update output
      # push!(ind_out,ind_vec_new);
      # push!(loc_out,loc_vec_new);
      push!(time_out,t);
      push!(N_out,N);

      #ERRORS
      #Break loop if extinction occurs
      if length(ind_vec_new) < 1
        println("Extinction has occured at t=",round(t,2)," and loop ",tic)
        break
      end
      #Break loop if resources go extinct and the other populations run away
      if NR == 0 && (NF + NH) > size
        println("Runaway growth has occured at t=",round(t,2)," and loop ",tic)
        break
      end

    end #end if Rate > 0

  end #end while loop over t

  println("Simulation successful at t= ",round(t,2)," and loop= ",tic)

  # return ind_out,loc_out,time_out,prop_out,N_out
  return time_out,prop_out,N_out
end #end function
