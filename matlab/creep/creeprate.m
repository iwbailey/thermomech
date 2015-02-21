function strainRate = creeprate( stress, A, n, E, T )
%CREEPSTRENGTH Get the strain rate at a given stress
%
% strainRate = creeprate( stress, A, n, E, T )

% Gas constant
R_g = 8.3144621;

% Temperature term
temperatureTerm = A.*exp( -E./(R_g * T));

% Go from stress to strain rate
strainRate = (stress.^n).*temperatureTerm;

end
