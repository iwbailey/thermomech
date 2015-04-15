function [stiff, selfstiff] = stiffnessmatrix( dx, dz, z0, mu, nl, nd )
%STIFFNESSMATRIX Calculate the slip to stress transformation for the fault
%
% stiff = stiffnessmatrix( dx, dz, z0, mu, nl, nd )
%
% dx along-strike width of a slip cell
% dz vertical dimension of a slip cell
% z0 depth at the center of the top cell
% mu rigidity
% nl = number of cells along strike
% nd = number of cells down

% Depth of cells 
z = z0 + ( (0:nd-1) )*dz;

% Difference in strike index btw source and receiver
diffx = (0:nl-1) * dx;

% Initialize the matrix
stiff = cell(1,nd);
selfstiff = zeros( nd, nl);

% Loop through all receiver depths
for i=1:nd,
    % Calculate strain contribution 
    stiff{i} = 2*mu*strikeslipstrain( 0.5*dz, 0.5*dx, z(i), diffx, z );
    selfstiff(i, :) = stiff{i}(i,1)*ones(1,nl);
end


end
