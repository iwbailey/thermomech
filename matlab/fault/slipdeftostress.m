function stress = slipdeftostress( slipdef, stiffnessmatrix, nx, nz )
%SLIPDEFTOSTRESS calculate stress on fault based on slip deficit
%
% stress = slipdeftostress( slipdef, stiffnessmatrix )

% dimensions
if( nargin < 4 ),
    [nz, nx] = size( slipdef );
end

% Loop through all cells
tic
stress = zeros(nz,nx);
for ik=1:nx,
    % along strike indices of stiffness matrix
    idx = [nx-(nx-ik):-1:2,  1:(nx-ik+1)];

    % Loop through receiver depths
    for j=1:nz,
        
        % Multiply slip def by stiffness matrix for this receiver cell
        stressFromCells = -slipdef.*stiffnessmatrix{j}(:,idx);
        
        % Sum to get stress
        stress( j, ik) = sum( stressFromCells(:) );
    end
end


end
