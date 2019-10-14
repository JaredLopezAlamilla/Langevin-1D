function [Y] = Murayama1D(diffu,ga,dt,dVx,X,wn)

  Y=X-dt/ga*dVx+sqrt(2*diffu)*wn;

end
