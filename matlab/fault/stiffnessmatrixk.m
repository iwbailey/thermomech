function [stiffk, selfstiff] = stiffnessmatrixk( dx, dz, z0, mu, nl, nd )
%STIFFNESSMATRIXK Calculate the slip to stress transformation for the fault in
% wavenumber domain
%
% stiff = stiffnessmatrixk( dx, dz, z0, mu, nl, nd )
%
% dx along-strike width of a slip cell
% dz vertical dimension of a slip cell
% z0 depth at the center of the top cell
% mu rigidity
% nl = number of cells along strike
% nd = number of cells down
%
% stiff = cell array of nd matrices. Each matrix is nd x nl

% Depth of cells 
z = z0 + ( (0:nd-1) )*dz;

% Difference in strike index btw source and receiver
diffx = (0:nl-1) * dx;

stiffM = zeros( nd, nl, nd);

% Loop through all receiver depths
for iDepth=1:nd,
    % Calculate strain contribution 
    stiffM(:,:,iDepth) = 2*mu*strikeslipstrain( 0.5*dz, 0.5*dx, z(iDepth), diffx, z );
end

% Initialize the cell array
stiffk = cell(1,nd);
selfstiff = zeros(nd, nl);
for iDepth=1:nd,
    % Selfstiffness
    selfstiff( iDepth,: ) = stiffM(iDepth, 1, iDepth)*ones(1,nl);
    
    % Build the filter for a stress calculation at cell iDepth, 1 due to slip
    % deficit at all cells on the fault
    stiff = squeeze(stiffM(:,:,iDepth));
    
    % Make a filter
    K = [fliplr(stiff(:,2:nl)), stiff];
    
    % Save the fourier transform
    stiffk{iDepth} = fft2(K);
end


end
