% Script generates creep parameters, generates files for input into programs and
% makes check plots

%% Clean up
clear
close all

%% Parameters

% Parameters for BZ1996
nx = 128; % Number of slip cells along strike
nz = 32; % Number of slip cells down dip
fLength = 70.0;
fWidth = 17.5;
zBD = 10; % Brittle-ductile transition depth in km
vPl = 35; % Plate velocity in mm/yr
cohesion = 60;
fs = 0.75; % Friction coefficient
dSigmaEff_dz = 180.0; % Change in effective normal stress (bar) per km
tauratioz = 4; % Factor decrease in creep strength from zBD to fault base

% Parameters for Arrhenius setup
actEnergy = 130e3; % Activation energy
Tsurface = 273+20; % Surface temperature

%dTdz = 30;
dTdz = 25;
R_g = 8.3144621;

ofileBZ = 'creepcoef.bz1996.in';
ofileArrh_creep = 'creepcoef.arrh.in';
ofileArrh_activ = 'activEnergy.arrh.in';

%% Get creep parameters on the fault using BZ1996 model
crpBZ = creep_bz1996( nx, nz, fLength, fWidth, zBD, tauratioz, vPl, ...
                      fs, dSigmaEff_dz );

%% Get new parameters
% Stress at zBD
stressBD = fs*(cohesion + dSigmaEff_dz*zBD);

% temperature at zBD
tempBD = Tsurface + dTdz*zBD;

[crpArrh, faultE] = creep_arrh( nx, nz, fLength, fWidth, actEnergy, zBD, ...
                                stressBD, tempBD, vPl, R_g, dTdz);

%% Write to file(s)
fid = fopen( ofileBZ, 'w');
fprintf(fid, '%14.10e\n', crpBZ);
fclose(fid);
fprintf('Written to %s\n', ofileBZ)

fid = fopen( ofileArrh_creep, 'w');
fprintf(fid, '%14.10e\n', crpArrh);
fclose(fid);
fprintf('Written to %s\n', ofileArrh_creep)

fid = fopen( ofileArrh_activ, 'w');
fprintf(fid, '%14.10e\n', faultE);
fclose(fid);
fprintf('Written to %s\n', ofileArrh_activ)

%% Calculate strength
xcell = (0.5 + (0:nx-1))*fLength/nx;
zcell = (0.5 + (0:nz-1)')*fWidth/nz;
faultTemp = repmat(Tsurface + zcell*dTdz, 1, nx);
brittleStrength = repmat( cohesion + dSigmaEff_dz*zcell, 1, nx );

% Vpl = crp * tau^3 * exp( -E/R*T )
% tau = (( Vpl/crp )*exp( E/R*T ) )
crpStrengthBZ = ((vPl./crpBZ).*exp( zeros(nz,nx)./(R_g*faultTemp) ) ).^(1/3);
strengthBZ = min( brittleStrength, crpStrengthBZ);

crpStrengthArrh = ((vPl./crpArrh).*exp( faultE./(R_g*faultTemp) ) ).^(1/3);
strengthArrh = min( brittleStrength, crpStrengthArrh);

%% Generate plots

% Plot creep coefficients
figure
imagesc(xcell, zcell, crpBZ);
cb = colorbar;
ylabel(cb, 'creep coefficient')
title('BZ 1996 creep coefficients')

figure
imagesc(xcell, zcell, crpArrh.*exp(-faultE./(R_g*faultTemp)) );
cb = colorbar;
ylabel(cb, 'creep coefficient')
title('Arrh creep coefficients')

% Plot strength
figure
imagesc( xcell, zcell, strengthBZ);
cb = colorbar;
ylabel(cb, 'strength [bar]')
title('BZ 1996')

% Plot strength
figure
imagesc( xcell, zcell, strengthArrh);
cb = colorbar;
ylabel(cb, 'strength [bar]')
title('Arrhenius Creep')

% Plot strength difference
figure
imagesc( xcell, zcell, strengthArrh-strengthBZ);
cb = colorbar;
ylabel(cb, '\Delta strength [bar]')
title('Arrhenius - BZ1996 Creep')


% Plot strength profile
icol = round(0.5*nx);
figure
plot( strengthBZ(:,63), zcell, '+k', 'linewidth', 1.5);
hold on
plot( strengthBZ(:,64), zcell, '+k', 'linewidth', 1.5);
plot( strengthBZ(:,63), zcell, '+k', 'linewidth', 1.5);

plot( strengthArrh(:,63), zcell, 'xb', 'linewidth', 1.5);
plot( strengthArrh(:,64), zcell, 'xb', 'linewidth', 1.5);
plot( strengthArrh(:,63), zcell, 'xb', 'linewidth', 1.5);
set( gca, 'ydir', 'reverse');
grid on;
xlabel('Strength [bar]')
ylabel('Depth [km]')