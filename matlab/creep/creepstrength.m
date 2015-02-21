function strength = creepstrength( strainRate, A, n, E, T )
%CREEPSTRENGTH Get the stress necessary for a given strain rate
%
% strength = creepstrength( strainRate, A, n, E, T )

% Gas constant
R_g = 8.3144621;

% Temperature term
temperatureTerm = A.*exp( -E./(R_g * T));

% Go from strain rate to stress
strength = (strainRate ./ temperatureTerm ).^(1/n);
end
