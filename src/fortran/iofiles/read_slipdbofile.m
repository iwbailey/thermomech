% Matlab/octave script for reading the output file

% Clean up variables
clear

%% Define parameters for script
ifile = 'slipdb2_test.unif.dtau12pm6'
nl = 128;
nd = 32;
dx = 70/nl;

%% Open file and read into variables
fid = fopen( ifile, 'rb');

tmp = fread( fid, 1, 'integer*4')

upl = fread( fid, 1, 'float64')/1e3
u1 = reshape( fread( fid, nl*nd, 'float64'), nl, nd)';

isOk = true;
iline = 1;
while(isOk)
    try
        tmp = fread( fid, 1, 'int64');

        nhypo(end+1) = fread( fid, 1, 'integer*4');
        ihypo(end+1) = fread( fid, 1, 'integer*4');
        jhypo(end+1) = fread( fid, 1, 'integer*4');
        thypo(end+1) = fread( fid, 1, 'real*8');
        nSlip(end+1) = fread( fid, 1, 'integer*4');
        avgSlip(end+1) = fread( fid, 1, 'real*8');

    catch ME
        fprintf('...%i eqks\n',numel(nhypo) );
        isOk = false;
        break
    end
end
fclose(fid);

p0 = nSlip*(1e3*dx)^2 .* (1e-3*avgSlip); % in m^3
m0 = 30e9*p0; % in Nm
mag = (log10(m0) - 9.1)./1.5;

%% Plot hypocenters
figure
plot( dx*(ihypo-0.5), dx*(jhypo-0.5), '.' );
set( gca, 'ydir' , 'reverse')
axis equal
ylim([0,nd*dx])
xlim([0,nl*dx])
title('Hypocenters')

%% Plot depth distribution
n = hist( jhypo, 1:nd );
n = [n(:)'; n(:)'];
depths = [ (0:nd-1)*dx; (1:nd)*dx];
figure
plot( n(:), depths(:), '-k', 'linewidth', 2);
xlabel('Count')
ylabel('Depth')
grid on
set( gca, 'ydir', 'reverse')

%% Plot time vs logM0
figure
plot( thypo, m0, '.k');
set( gca, 'yscale', 'log');
grid on
xlabel('Year')
ylabel('M0')

xlim([150,300])

%% Plot magnitude frequency distribution
figure
plot( sort(mag,'descend'), 1:numel(p0), '-ok', 'linewidth', 2 );
set( gca, 'yscale', 'log');
grid on
xlabel('M_W')
ylabel('Number(m\geq M_W)');

%% Plot distribution of interevent times
figure
hold on
colors = {'k', 'b', 'g', 'r'};
minm = 3:6;
for i=1:numel(minm),
    plot( sort( diff(thypo(mag>=minm(i))), 'descend'), 1:nnz(mag>=minm(i))-1, '-o', ...
          'linewidth', 2, 'markersize', 3, 'color', colors{i});
end
set( gca, 'xscale', 'log', 'yscale', 'log');
grid on
xlabel('Inter-event time')
ylabel('Number');
legend('M_W\geq 3', 'M_W\geq 4', 'M_W\geq 5', 'M_W\geq 6')
