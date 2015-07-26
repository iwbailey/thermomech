% Matlab/octave script for reading the output file

% Clean up variables
%clear

%% Define parameters for script
ifile = 'lastdb2.unif.t150.dtau12pm6'
nl = 128;
nd = 32;
dx = 70/nl;

%% Open file and read into variables
fid = fopen( ifile, 'rb');

tmp = fread( fid, 1, 'integer*4')

upl = fread( fid, 1, 'float64')/1e3
u1 = reshape( fread( fid, nl*nd, 'float64'), nl, nd)';

fclose(fid);

%% Plot u1
figure
imagesc( dx*((1:nl)-0.5), dx*((1:nd)-0.5), u1 )
cb = colorbar;
ylabel( cb, 'u1 [mm]')
title('Slip after initialization program')
return
