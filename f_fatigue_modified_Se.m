function [ Se, ka, kb, kc ] = f_fatigue_modified_Se(surface,d_eq,load_type,Sut,Sepr)
% f_fatigue_modified_Se : calculate Se = ka·kb·kc·Se'
%   ka: surface factor (options: ground, machined, hotrolled, forged, casted)
%   kb: size factor (introduce equivalent diameter d_eq in 'mm')
%   kc: load factor (options: bending, axial, torsion)
%       NOTE: for 'torsion' kc = 1 is only valid when the fatigue analysis
%       is done with the von Mises equivalent stress. For pure torsion
%       analysed with shear stresses, Shigley gives kc = 0.59.
%   Sut: ultimate tensile strength in 'MPa'
%   Sepr: Se', introduce Se' value in 'MPa' or material for estimation
%         options: steel, castiron, aluminium, copper, polymer

% --- Surface factor ka = a·Sut^b ---
switch surface
    case 'ground'
        a=1.58; b=-0.085;
    case 'machined'
        a=4.51; b=-0.265;
    case 'hotrolled'
        a=57.7; b=-0.718;
    case 'forged'
        a=272;  b=-0.995;
    case 'casted'
        a=1;    b=0;
    otherwise
        error(['f_fatigue_modified_Se: surface option not valid: ''%s''. ', ...
               'Options: ground, machined, hotrolled, forged, casted'], surface);
end
ka=a*Sut^b;

% --- Size factor kb ---
switch load_type
    case 'axial'
        kb=1;
    case {'bending','torsion'}
        if (d_eq<2.79)
            kb=(2.79/7.62)^(-0.107);
        elseif (d_eq>51)
            kb=1.51*d_eq^(-0.157);
        else
            kb=(d_eq/7.62)^(-0.107);
        end
    otherwise
        error(['f_fatigue_modified_Se: load_type option not valid: ''%s''. ', ...
               'Options: bending, axial, torsion'], load_type);
end

% --- Load factor kc ---
switch load_type
    case 'bending'
        kc=1;
    case 'axial'
        kc=0.85;
    case 'torsion'
        kc=1;      % only valid with von Mises equivalent stress (see NOTE)
end

% --- Endurance limit of the test specimen Se' ---
if ischar(Sepr)
    switch Sepr
        case 'steel'
            if (Sut<1400) Sepr=0.5*Sut; else Sepr=700; end
        case 'castiron'
            if (Sut<400) Sepr=0.4*Sut; else Sepr=160; end
        case 'aluminium'
            if (Sut<330) Sepr=0.4*Sut; else Sepr=130; end
        case 'copper'
            if (Sut<280) Sepr=0.4*Sut; else Sepr=100; end
        case 'polymer'
            Sepr=0.2*Sut; % to Sepr=0.4*Sut;
        otherwise
            error(['f_fatigue_modified_Se: Sepr option not valid: ''%s''. ', ...
                   'Introduce a value in MPa or one of: steel, castiron, ', ...
                   'aluminium, copper, polymer'], Sepr);
    end
end

Se=Sepr*ka*kb*kc;
end
