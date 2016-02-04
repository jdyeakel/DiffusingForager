function starvingforager_event_nodiff(L,dim,initsize,t_term,alpha,K,sigma,rho,lambda,mu)
  #Read in packages/function
  #ipbc :: torus movement
  #include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/ipbc.jl")


  #Initiate values
  #Lattice dimension

  S = (L-2)^dim;


  #Initial state values
  ind_vec = [rand(collect(1:2),Int(2*round(initsize/3)));zeros(Int64,Int(round(initsize/3)))];
  #OR Pure Resource Start (for testing)
  #ind_vec = zeros(Int64,initsize);


  #Initial location values
  #2/3 of initsize will be foragers; 1/3 will be resources
  #replace = false because want only one resource per site
  resourceloc_vec = sample(collect(1:S),Int(round(initsize/3)),replace=false);
  loc_vec = [rand(collect(1:S),Int(2*round(initsize/3)));resourceloc_vec];

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
  noresourcesites = collect(1:S);
  #Deletes locations of sites with resources from the list of open sites (noresourcesites)
  deleteat!(noresourcesites,sort(resourceloc_vec));

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
  F = NF/S;
  H = NH/S;
  R = NR/S;
  prop_out = Array{Float64}(3,1);
  prop_out[:,1] = [F,H,R];
  push!(N_out,N);
  #push!(prop_out,copy(prop));

  tic = 0;

  # F_pr_line = 0.0;
  # H_pr_line = 0.0;
  # R_pr_line = 0.0;

  # id = Array(Float64,1);

  while t < (t_term-1)
    tic = tic + 1;


    # #ERRORS
    # if length(ind_vec) < 1
    #   print("Extinction has occured at t=",round(t,2)," and loop ",tic)
    #   break
    # end

    #Calculate Rate
    Rate = F*(lambda + sigma*(K-R)) + H*(rho*R + mu) + R*(alpha*(K-R) + (F+H)); # (1 - (N/S)) +

    # TESTING
    # Rate = (NF/N)*(lambda + sigma*(K-R)) + (NH/N)*(rho*R + mu) + (NR/N)*(alpha*(K-R) + (F+H));

    dt = 1/(Rate*N);
    if Rate == 0
      println("Welcome to Daisy World")
      break
    end

    #Construct probability lines, which are a function of R, S, F
    #Grow <-----> Starve <-----> Diffuse

    F_pr_line = lambda/Rate;
    #F_pr_line[2] = copy(F_pr_line[1]) + (sigma*(K-R))/Rate;

    #Recover <-----> Mortality <-----> Diffuse

    H_pr_line = (rho*R)/Rate;
    #H_pr_line[2] = copy(H_pr_line[1]) + mu/Rate;

    #Grow <-----> Consumed
    R_pr_line = (alpha*(K-R))/Rate;
    # R_pr_line[2] = R_pr_line[1] + (F+H)/Rate;



    #Randomly select an individual (R,S,F) with probability 1/N
    #ind thus represents the POSITION of the individual
    #Update total number of individuals

    #Update the total

    #Choose a random individual position
    id = rand(collect(1:N));

    #Also, we need to make an ind_vec_old that is NOT updated for the following sequence of IF statements
    ind_vec_old = copy(ind_vec);

    #If randomly selected individual is Full
    if ind_vec_old[id] == 2
      state = 2;

      draw_event = rand();

      #GROW
      if draw_event < F_pr_line
        push!(ind_vec,2); #ensure the new individual is in the full state
        #Offspring appears at a random location
        location = rand(collect(1:S));
        push!(loc_vec,location);
        NF = NF+1;
      end

      #STARVE
      if draw_event >= F_pr_line #&& draw_event < F_pr_line[2]
        ind_vec[id] = 1;

        NH = NH+1;
        NF = NF-1;
      end

      # #MOVE
      # if draw_event >= F_pr_line[2]
      #   #Move to a random location
      #   loc_vec[id] = rand(collect(1:S));
      # end


    end



    #If randomly selected individual is Hungry
    if ind_vec_old[id] == 1
      state = 1;

      #Draw a random event
      #Recover, die, or move??
      draw_event = rand();

      #Recover (become full)
      if draw_event < H_pr_line
        ind_vec[id] = 2;
        NH = NH-1;
        NF = NF+1;
      end

      #Die
      if draw_event >= H_pr_line #&& draw_event < H_pr_line[2]
        deleteat!(ind_vec,id);
        deleteat!(loc_vec,id);
        NH = NH - 1;
      end

      # #MOVE
      # if draw_event >= H_pr_line[2]
      #   #Move to a random location
      #   loc_vec[id] = rand(collect(1:S));
      # end

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

        #Choose random position on the noresourcesite list
        newresourcepos = rand(collect(1:length(noresourcesites)));
        #Define the new position as the empty site for resource growth
        location = noresourcesites[newresourcepos];

        #Append the resource's location to the END of the vector
        push!(loc_vec,location);
        #Update Tally
        NR = NR + 1;

        #Update the no resource site by deleting the position that has been filled
        deleteat!(noresourcesites,newresourcepos);

      end

      #BECOME CONSUMED!
      if draw_event >= R_pr_line

        #Update the no resource site by adding the position that is now empty
        #Needs to be done BEFORE the ind/loc information is deleted.
        push!(noresourcesites,copy(loc_vec[id]));

        #Delete individual from the list
        deleteat!(ind_vec,id);
        #Delete individual location from the list
        deleteat!(loc_vec,id);
        NR = NR - 1;

      end

    end

    #Recalculate the total # of the foragers and resources
    N = NF + NH + NR;
    #Recalculate the densities of each
    F = NF/S;
    H = NH/S;
    R = NR/S;

    prop = [F,H,R];

    #push!(prop_out,copy(prop));
    prop_out = hcat(prop_out,prop);

    #Advance time
    t = t + dt;

    if t > next_time_int
      println("time= ",round(t,2))
      next_time_int = round(t+1,0)
    end

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
    if NR == 0 && (NF + NH) > S
      println("Runaway growth has occured at t=",round(t,2)," and loop ",tic)
      break
    end



  end #end while loop over t

  println("Simulation successful at t= ",round(t,2)," and loop= ",tic)

  # return ind_out,loc_out,time_out,prop_out,N_out
  return time_out,prop_out,N_out
end #end function
