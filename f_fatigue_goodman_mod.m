function [ SF_f ] = f_fatigue_goodman_mod(sigma_m,sigma_a,Se,Sut)
% f_fatigue_goodman_mod : fatigue safety factor with the modified Goodman
%   criterion, assuming a proportional load line (sigma_a/sigma_m constant).
%
%   sigma_m : mean stress in 'MPa' (with sign: < 0 compression)
%   sigma_a : alternating stress in 'MPa' (>= 0)
%   Se      : modified endurance limit in 'MPa' (e.g. f_fatigue_modified_Se)
%   Sut     : ultimate tensile strength in 'MPa'
%
%   SF_f : fatigue safety factor
%          sigma_m >= 0 : SF_f = 1 / (sigma_a/Se + sigma_m/Sut)
%          sigma_m <  0 : SF_f = Se / sigma_a  (compressive mean stress is
%                         not taken into account)

if sigma_a < 0
    error('f_fatigue_goodman_mod: sigma_a must be >= 0 (sigma_a = %g)', sigma_a);
end
if sigma_a == 0 && sigma_m <= 0
    error('f_fatigue_goodman_mod: sigma_a = 0 and sigma_m <= 0, there is no fatigue load');
end

if sigma_m >= 0
    SF_f = 1 / (sigma_a/Se + sigma_m/Sut);
else
    SF_f = Se / sigma_a;
end
end
