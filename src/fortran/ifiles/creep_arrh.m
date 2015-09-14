function [faultCrp, faultE] = creep_arrh( nx, nz, fLength, fWidth, actEnergy, ...
                                          zBD, stressBD, tempBD, vPl, R_g, dTdz)
% Setup the creep parameters for the arrhenius verion of creep

% Distance from BD transition to use for perturbation
dzPerturb = 1.25;

% Width of border region that is not adjusted
borderWidth = 3.75;

% Probability of adjustment
probAdjust = 0.2;

% Get slip cell dimensions
dx = fLength/nx; % Slip cell length in km
dz = fWidth/nz; % Slip cell depth in km
xBD = fWidth-zBD;

% Arrhenius constant/fault width
%vPl = faultWidth * const * sigma.^3 * exp( -actEnergy / (R_g*tempProfile) )
%vPl = const * sigma.^3 * exp( -actEnergy / (R_g*tempProfile) )
crpConst = (vPl/stressBD^3) * exp( actEnergy / (R_g*tempBD) )

% Expand to fault
faultE = actEnergy*ones(nz,nx);
faultCrp = crpConst*ones(nz,nx);

% Add temp-independent creeping boundaries
[faultCrp, faultE] = faultcrpboundaries( faultCrp, faultE, dx, dz, zBD, ...
                                         tempBD, dTdz, R_g );

% Add the randomness
[faultCrp, isRand] = add_randomness( faultCrp, dx, dz, crpConst, actEnergy, ...
                                     zBD, tempBD, dTdz, dzPerturb, borderWidth, ...
                                     probAdjust, R_g);

% Correction for the boundary regions
faultTemp = get_faulttemperature( dTdz, tempBD, zBD, dz, nz, nx );
isAdjust = isRand & faultCrp ~= crpConst & faultE == 0;
faultCrp(isAdjust) = faultCrp(isAdjust).*...
    exp( -actEnergy ./ (R_g*faultTemp(isAdjust) ) );

end

%------------------------------------------------------------------------------
function [crp, actEnergy] = faultcrpboundaries( crp, actEnergy, dx, dz, ...
                                                zBD, tempBD, dTdz, R_g)
% Apply creeping boundaries of the bottom to the sides
[nz, nx] = size(crp);

xBD = nz*dz - zBD;

% Define cell positions
xcell = repmat( dx*(0.5+(0:nx-1)) , nz, 1);
zcell = repmat( dz*(0.5+(0:nz-1)'), 1, nx);

% Get current creep profile on its side
% Tprofile = T0 + z.dTdz
%          = (Tbd - zBD.dTdz) + z.dTdz
%          = Tbd + dTdz*(z-zBD)
xBDfar = nx*dx - xBD;
sideDistFar = max( xcell, nx*dx-xcell);
Txsect = tempBD + dTdz*(sideDistFar - xBDfar);
crpXsect = crp.*exp( -actEnergy./(R_g*Txsect) );


% Define boundaries
sideDist = min( xcell, nx*dx-xcell);
bottomDist = nz*dz-zcell;
isBoundary = sideDist < xBD & sideDist < bottomDist;

% Replace parameters
crp(isBoundary) = crpXsect(isBoundary);
actEnergy(isBoundary) = 0.0;

end

%------------------------------------------------------------------------------
function [crp, isHit] = add_randomness( crp, dx, dz, crpConst, actEnergy, zBD, tempBD, ...
                                        dTdz, dzPerturb, borderWidth, probAdjust, R_g)

% Fault dimensions
[nz, nx] = size(crp);

% Get temp on the fault
zcell = repmat( dz*(0.5+(0:nz-1)'), 1, nx);

% Get creep parameter at brittle and creep depths
Tbrittle = tempBD + dTdz*(-dzPerturb);
Tductile = tempBD + dTdz*(dzPerturb);
%crpBrittle = crp * exp( -actEnergy / (R_g*Tbrittle ) );
%crpDuctile = crp * exp( -actEnergy / (R_g*Tduct ) );

% Set the random seed
rand( "state", 1.2345*(1:625)' );

% Fault dimensions
faultDepth = nz*dz;
faultLength = nx*dx;

% Calculate depths at the center of each cell
zcell = repmat( 0.5*dz + dz*(0:nz-1)', 1, nx);

% Horizontal positions at center of each cell
xcell = repmat( 0.5*dx + dx*(0:nx-1), nz, 1);

% Temperature on fault
faultTemp = tempBD + dTdz*(zcell - zBD);

% Generate random numbers
r = rand( nz, nx);

% Keep 20% of the values within bounds
isHit = r <= probAdjust & ...
        zcell < (faultDepth - borderWidth) & ...
        xcell >= borderWidth & xcell < (faultLength - borderWidth);

% Identify creeping parts of the fault
isCreep = zcell >= zBD;

% Give creeping parts a brittle strength
crp( isHit & isCreep ) = crpConst .* ...
    exp( actEnergy * ( 1./(R_g*faultTemp(isHit&isCreep)) - 1./(R_g*Tbrittle) ...
                       ) );

% Give brittle parts a creep strength
crp( isHit & ~isCreep) = crpConst .* ...
    exp( actEnergy * ( 1./(R_g*faultTemp(isHit&~isCreep)) - 1./(R_g*Tductile) ...
                       ) );

end

%-------------------------------------------------------------------------------
function faultTemp = get_faulttemperature( dTdz, tempBD, zBD, dz, nz, nx );

zcell = repmat( dz*(0.5+(0:nz-1)'), 1, nx);
faultTemp = tempBD + dTdz*(zcell - zBD);

end
