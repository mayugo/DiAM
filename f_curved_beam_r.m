function [rc,rn,A,I]=f_curved_beam_r(ri,ro,b)
% Variables d'entrada (segons figura): ri, ro, b (amplada)
% Variables de sortida: radi centroide (rc), radi fibra neutre (rn), 
%                       àrea (A), i moment de segon ordre (I)
h = ro-ri;  
A = h*b;

rc= ri + (h/2);
rn= h/(log(ro/ri));

I = b*h^3/12  ;
end