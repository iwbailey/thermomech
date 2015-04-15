function m0 = scalarmoment( slipdef1, slipdef2, cellArea, rigidity )
%SCALARMOMENT calculate scalar moment based on slip deficit change
%
% m0 = scalarmoment( slipdef1, slipdef2, cellArea, rigidity )

m0 = rigidity*scalarpotency( slipdef1, slipdef2, cellArea );


end