function creepCoeffs = faultcreep_bz96( nL, nD, cellLength, cellHeight, xBD, zBD, ...
                                        strainRate, strength_zBD )
%FAULTCREEP_BZ96 Generate A values that are like the Ben-Zion 1996 Creep mask
%
% faultA = faultcreep_bz96( nL, nD, cellL, cellH, xBD, zBD, strainRate, strength_ZBD )
%

% Compute the distance from the side  for all cells
if( mod( nL, 2) == 0 ),
    sideDist = cellLength*([1:0.5*nL, 0.5*nL:-1:1]-0.5);
else
    sideDist = cellLength*([1:0.5*(nL+1), 0.5*(nL-1):-1:1]-0.5);
end
sideDist = repmat( sideDist, nD, 1);

% Compute distance from the bottom
bottomDist = repmat( cellHeight*((nD:-1:1)'-0.5), 1, nL);

% Keep the side or bottom result, whichever is lowest strength and highest
% creep coeff
creepCoeffs = max( ...
    creep_bz96( sideDist, strainRate, strength_zBD, xBD ), ...
    creep_bz96( bottomDist, strainRate, strength_zBD, zBD) );

% TODO: assign brittle strength to 20% of creep cells

end

% ------------------------------------------------------------------------------
function c = creep_bz96( x, vPl, strengthz, xBD )
%*
% Return the temperature independent creep parameter in Ben-Zion
% (1996), JGR, eqn 2

% x = distance from edge of fault (note difference from BZ1996 where it is along strike coord)
% vPl = plate velocity
% strengthz = background strength at that depth
% xBD = creeping boundary width


n = 3.0; % stress exponent

% Stress ratio of 4 means stress required for vPl at x=0 is 0.25 strength at xBD
ln_tauRatio = log(4.0);

c = ( vPl/strengthz^n ) * exp( n*ln_tauRatio*(1 - x/xBD) );

end