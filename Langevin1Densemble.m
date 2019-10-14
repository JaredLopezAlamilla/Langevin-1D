%-------------------------------------------------------------------------
%    This code reads the computes an ensemble of single-molecule trajectory
%      of overdamped Brownian motion on a fre-energy landscape
%               using the overdamped Langevin equation
%-------------------------------------------------------------------------

% ---- Thermal bath Parameters ----
kT=1; diffu=1; ga=kT/diffu; L=1;
% ---- Simulation Parameters ----
 dt=0.4e-2; Time=650; Ensemble=240;
NT=ceil(Time/dt);
rng('shuffle'); 
wn=sqrt(dt)*normrnd(0,1,[Ensemble,NT]);
% ---- Potential parameters ----
A=2;
% ---- We need to declare the potential as a symbolic function 1st
syms x; fx=2;
V=A*cos(2*pi*x/L)-fx*x;
dVx=diff(V,x);
FFx=matlabFunction(dVx); 
VV=matlabFunction(V);
% ---- Initialization of simulation ----
XX=zeros(Ensemble,NT+1);
XX(:,1)=normrnd(0,1,[Ensemble,1]);  %--- random initial position in x

for l=1:NT
    [Y] = Murayama1Densemble(diffu,ga,dt,FFx(XX(:,l)),XX(:,l),wn(:,l),Ensemble);
    XX(:,l+1)=Y(:,1);
end

TT=0:dt:Time;
if length(TT)<length(XX)
    TT=[TT TT(end)+dt];
else
end
% ---- Properties of simulated Brownian motion ----
Ppos=zeros(Ensemble,100);
[Ppos(1,:),Xedges,Yedges] = histcounts(mod(XX(1,:),2*pi),100);

for kk=2:Ensemble
    [Ppos(kk,:),~,~] = histcounts(mod(XX(kk,:),2*pi),100);
end

Ppos=mean(Ppos,1);

centers=(circshift(Xedges',1)-Xedges')/2+Xedges';
centers=centers(2:end,1)';
Ppos=Ppos/trapz(centers,Ppos); % Steady-state probability

P2 = polyfit(TT,mean(XX(:,:)-XX(:,1)),1);
vx=P2(1); % 1st moment: drift velocity


DDeff=mean(var(XX).')'./(2*dt)/NT; % 2nd moment: diffusion






