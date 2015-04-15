function p0 = scalarpotency( slipdef1, slipdef2, cellArea )
%SCALARPOTENCY calculate scalar potency based on slip deficit change
%
% p0 = scalarpotency( slipdef1, slipdef2, cellArea )

slip = slipdef1 - slipdef2;
p0 = sum(slip(:)).*cellArea;


end