function [N,mode] = f_fatiga_limitada(sigma_m, sigma_a,S_ut,S_0,N_0,S_e,N_e)
% Calcula la vida esperada per a un estat de carrega equivalent segons la
% funcio de Basquin: Seqv = a N^b
%   sigma_m = tensio equivalent mitjana (nomes per sigma_m positives)
%   sigma_a = tensio equivalent alterna
%   S_ut    = resistencia a traccio
%   S_0     = resistencia a N_0, per R=-1
%   N_0     = cicles
%   S_e     = limit de fatiga a N_e, per R=-1
%   N_e     = cicles

if (sigma_m<0) % comprova que tensio mitjana es positiva
    N=0;
    mode = 'error: tensio mitjana negativa';
    return
end

sigma_max = sigma_m + sigma_a;
A = sigma_a/sigma_m;    % A= sigma_a  /sigma_m = 1-R / 1+R
if A==Inf
    Sfe=S_e;
    Sf0= S_0;
else
    Sfe = (1+A)*S_ut/(1+A*S_ut/S_e);
    Sf0 = (1+A)*S_ut/(1+A*S_ut/S_0);
end

if sigma_max>S_ut
    N=0;
    mode = 'falla estatica';
elseif sigma_max>Sf0
    b_1=(log10(S_ut)-log10(Sf0))/(log10(1)-log10(N_0));
    a_1=10^(log10(S_ut)-b_1*log10(1));
    N=(sigma_max/a_1)^(1/b_1);
    mode = 'baix #cicles: vida limitada';
elseif sigma_max>Sfe
    b_2=(log10(Sf0)-log10(Sfe))/(log10(N_0)-log10(N_e));
    a_2=10^(log10(Sf0)-b_2*log10(N_0));
    N=(sigma_max/a_2)^(1/b_2);
    mode = 'alt #cicles: vida limitada';
else
    N=Inf;
    mode = 'alt #cicles: vida ilimitada';
end
end
