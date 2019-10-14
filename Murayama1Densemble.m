function [Y] = Murayama1Densemble(diffu,ga,dt,dVx,X,wn,Ensemble)

	Iden=spdiags(ones(Ensemble,1),0,Ensemble,Ensemble);
    
	diffuX=spdiags(sqrt(2*diffu)*ones(Ensemble,1),0,Ensemble,Ensemble);
    gaX=spdiags((-dt/ga)*ones(Ensemble,1),0,Ensemble,Ensemble); 
        
    Y=Iden*X+gaX*dVx+diffuX*wn;
    
end
