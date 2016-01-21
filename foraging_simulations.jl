#include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event.jl")
include("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event.jl")
using StatsBase
using Gadfly


L = 50;
dim = 2;
initsize = (L-2)^dim;
t_term = 100;
alpha = 0.5;
K = 1;
sigma = 0.6;
rho = 0.25;
lambda = 0.2;
mu = 0.2;
DF = 0.2;
DH = 0.2;

ind_out, loc_out, time_out, prop_out = starvingforager_event(L,dim,initsize,t_term,alpha,K,sigma,rho,lambda,mu,DF,DH);



#Posthoc analysis
steps = length(ind_out);
pop_F = zeros(steps);
pop_H = zeros(steps);
pop_R = zeros(steps);
for i=1:steps
  pop_F[i] = length(find(x->x==2,ind_out[i]));
  pop_H[i] = length(find(x->x==1,ind_out[i]));
  pop_R[i] = length(find(x->x==0,ind_out[i]));
end



plot(layer(x=time_out,y=pop_F,Geom.line,Theme(default_color=color("green"))),
layer(x=time_out,y=pop_H,Geom.line,Theme(default_color=color("orange"))),
layer(x=time_out,y=pop_R,Geom.line,Theme(default_color=color("blue"))))


plot(x=time_out,y=prop_out)
