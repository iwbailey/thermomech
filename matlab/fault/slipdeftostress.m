function stress = slipdeftostress( slipdef, stiffnessmatrix, nx, nz )
%SLIPDEFTOSTRESS calculate stress on fault based on slip deficit
%
% stress = slipdeftostress( slipdef, stiffnessmatrix )

% dimensions
if( nargin < 4 ),
    [nz, nx] = size( slipdef );
end

% Loop through all cells
stress = zeros(nz,nx);
for ik=1:nx,
    % along strike indices of stiffness matrix
    idx = [nx-(nx-ik):-1:2,  1:(nx-ik+1)];
    for j=1:nz,
        stress(j,ik) = sum( flat( -slipdef.*squeeze(stiffnessmatrix(j,:,idx)) ) );
    end
end

end
%-------------------------------------------------------------------------------
function x = flat(x)
%FLAT flatten a matrix
x = x(:);
end