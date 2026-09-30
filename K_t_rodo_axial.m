function [ K_t] = K_t_rodo_axial( D_d , r_d, i_flag )
%Calcula la K_t d'un rodó sotmés a tracció amb un canvi de secció de D a d:
%   D = diàmetre gran; d= diàmetre petit; r=radi de l'entalla
%   D_d = D/d ; r_d = r/D
%   i_flag = 1, si es vol K_t per tensió normal
%   i_flag <> 1, si es vol K_t per von Mises

if i_flag==1 
    K_t = 0.493+0.480*D_d^(-2.43)+r_d^(-0.48)*...
      sqrt((3.43-3.41*D_d^2-0.0232*D_d^4)/(1-8.85*D_d^2-0.078*D_d^4));
else 
    K_t = 0.496+0.472*D_d^(-2.85)+r_d^(-0.48)*...
      sqrt((2.921-2.945*D_d^2-0.0217*D_d^4)/(1-9.59*D_d^2-0.053*D_d^4));
end
end