function stiff = stiffnessmatrix( dx, dz, z0, mu, nl, nd )
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

% Depth of cells */
z = z0 + ( (1:nd)+0.5 )*dz;

% Difference in strike index btw source and receiver
diffx = (1:nl) * dx;

% Initialize the matrix
stiff = zeros( nd, nd, nl);

for i=1:nd,
    stiff(i,:,:) = 2*mu*strikeslipstrain( 0.5*dz, 0.5*dx, z(i), ...
                                          diffx, z );
end


end
