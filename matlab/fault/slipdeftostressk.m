function stress2 = slipdeftostressk( slipdef, stiffnessmatrixk, nx, nz )
%SLIPDEFTOSTRESS calculate stress on fault based on slip deficit
%
% stress = slipdeftostress( slipdef, stiffnessmatrix )

% dimensions
if( nargin < 4 ),
    [nz, nx] = size( slipdef );
end

% FFt of slip def
slipdefk = fft2([flipud(-slipdef), zeros(nz,nx-1)]);

% Loop through all cells
stress2 = zeros(nz,nx);
for j=1:nz,
    tmp = ifft2( stiffnessmatrixk{j}.*slipdefk);
    stress2( j, :) = tmp(nz,nx:2*nx-1);
end

%stress2 = rot90(stress2,2);

end
