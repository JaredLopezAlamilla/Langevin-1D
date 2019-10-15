# Langevin 1D Simulations

This project contains the scripts neccesary to simulate 1D singe-molecule trajectories.
The dynamics of the system simulated is an overdamped Brownian motion over a periodic potential.
The sripts are written for matlab.

Functions required: Murayama2D.m, Langevin2D.m

Examples of use:
  
`>> Langevin2D;`

Plotting the results:

`>> figure(); plot(TT,XX(1,:));title('$x(t)$');`
