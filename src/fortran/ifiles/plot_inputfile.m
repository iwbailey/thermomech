% Use octave or matlab to plot the input files

% This script is basically just a reminder of how the arrays are ordered

ifile = 'creepcoef.bz1996.in';
nl = 128;
nd = 32;

% Read input file into a single column
t = load(ifile);

% reorder
t = reshape( t, nd, nl );

% Plot
imagesc( t );
cb = colorbar;
