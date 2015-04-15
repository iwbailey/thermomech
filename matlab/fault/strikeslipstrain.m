function strain = strikeslipstrain( halfwidth, halfdepth, zStrain, xStress, zStress )
%STRIKESLIPSTRAIN calc shear strain for point source due to unit slip
%
% strain = strikeslipstrain( xw2, xd2, x3, y1, y3 )
%
%     Chinnery solution for stress due to strike slip on a vertical, rectangular
%     fault patch  (following materials modified from an original coding by W.
%     D. Stuart provided Oct.'90 by T. E. Tullis)
%
% Copied from Ben-Zion Fortran code
%     Ref:
%     Chinnery, bssa, 1963
%
%
% xw2, xd2 are cell halfwidth and halfdepth
% tensor strain on fault plane at x1(=0),x3 (source cell center);
% due to stress at y1,y3 (stress calculation point)
%
% Multiple calculation points can be specified

% Get all coordinates for output
[xStress,zStress] = meshgrid(xStress, zStress);

% Follow summation/notation of paper
strain = ...
    + f(  halfwidth, zStrain+halfdepth, xStress, zStress ) ...
    - f(  halfwidth, zStrain-halfdepth, xStress, zStress ) ...
    - f( -halfwidth, zStrain+halfdepth, xStress, zStress ) ...
    + f( -halfwidth, zStrain-halfdepth, xStress, zStress );
end
%-------------------------------------------------------------------------------
function result = f( x1, x3, y1, y3)
% function does calculation

p = x3 + y3;
q = x3 - y3;
t = x1 - y1;

% Level 2 notation
s1 = sqrt( t.^2 + q.^2 );
s2 = sqrt( t.^2 + p.^2 );

% Level 3 notation
s1q = s1 + q;
s2p = s2 + p;

% Level 4 notation
f1 = 1./(s1.*s1q) + 1./(s2.*s2p);

f2 = (s2/4.0 + q)./(s2.*s2p.^2) - ...
    (p.^2 - q.^2).*(2*s2 + p)./(2*s2.^3 .*s2p.^2);

f4 = q./(s1.*(s1+t)) + p./(s2.*(s2+t));

% Combine to get result
result = ((2/3)*t.*(f1 + f2) + 0.5*f4)/(4*pi);
        
end