function mag = potmagnitude( scalarpotency )
%POTMAGNITUDE calculate the magnitude from the scalar potency
%
% Use the empirical relation of Ben-Zion and Zhu (GJI, 2002), quadratic form 
%
% m = potmagnitude( scalarpotency )

% Quadratic terms
a = 0.0612;
b = 0.988;

% Make sure we have the correct units for the scaling relation
km = 1e3;
cm = 1e-2;
c = -4.87 - log10( scalarpotency /(km*km*cm));

% Use the quadratic formula
mag = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);

end