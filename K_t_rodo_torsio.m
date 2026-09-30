function [ K_t] = K_t_rodo_torsio( D_d , r_d)
%Calcula la K_t d'un rodó sotmés a torsió amb un canvi de secció de D a d:
%   D = diàmetre gran; d= diàmetre petit; r=radi de l'entalla
%   D_d = D/d ; r_d = r/D
 
    K_t = 0.78+0.2*D_d^(-10)+r_d^(-0.46)*...
      sqrt((0.002-0.125*D_d^2+0.123*D_d^4)/(1-2.75*D_d^2+2.55*D_d^4));

end