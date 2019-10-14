%-------------------------------------------------------------------------
%    This code reads the computes an ensemble of single-molecule trajectory
%      of overdamped Brownian motion on a fre-energy landscape
%               using the overdamped Langevin equation
%-------------------------------------------------------------------------

% ---- Thermal bath Parameters ----
kT=1; diffu=1; ga=kT/diffu; L=1;
% ---- Simulation Parameters ----
 dt=0.4e-2; Time=650;
NT=ceil(Time/dt);
rng('shuffle'); 
wn=sqrt(dt)*normrnd(0,1,[1,NT]);
% ---- Potential parameters ----
A=2;
% ---- We need to declare the potential as a symbolic function 1st
syms x; fx=2;
V=A*cos(2*pi*x/L)-fx*x;
dVx=diff(V,x);
FFx=matlabFunction(dVx); 
VV=matlabFunction(V);
% ---- Initialization of simulation ----
XX=zeros(1,NT+1);
XX(:,1)=normrnd(0,1,[1,1]);  %--- random initial position in x

for l=1:NT
    [Y] = Murayama1D(diffu,ga,dt,FFx(XX(l)),XX(l),wn(l));
    XX(l+1)=Y(1);
end

TT=0:dt:Time;
if length(TT)<length(XX)
    TT=[TT TT(end)+dt];
else
end


