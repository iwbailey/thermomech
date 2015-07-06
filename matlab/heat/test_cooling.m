%
% test_cooling.cpp ---
%
% Output the cooling over a period of time after generating some heat via slip
% on a single fault patch
%
% Filename: test_cooling.cpp
% Description:
% Author: Iain W. Bailey
% Created: Wed Sep 25 11:31:04 2013 (-0400)
% Version: 1
% Last-Updated: Wed Sep 25 11:33:19 2013 (-0400)
%           By: Iain W. Bailey
%     Update #: 6
%
clear 
close all
Units = si_units();

%% Parameters
% PREM density of crustal rock kg/km^3
density = 2.6*Units.g/(Units.cm*Units.cm*Units.cm);

% Fault Width
faultWidth = 10.0*Units.cm;

% Diffusivity of rock, don't know origin 
diffusivity = 1.E-2*Units.cm*Units.cm/Units.s;

% Specific heat capacity for rock don't know origin
specHeat = 790*Units.J/(Units.kg*Units.K);

% Set the required variables
dx = (70.0*Units.km)/128; % width of a slip cell
dz = (17.5*Units.km)/32; % height of a slip cell
slip = 1*Units.m;  % amount of slip
stress = 100*Units.MPa; % shear stress during slip
creepSlipTime = 1*Units.hr;  % time for creep slip
creepSlipTime2 = 1*Units.day; % time for second creep slip
vRupture = 3.0*Units.km/Units.s; % rupture velocity
eqkSlipTime = slip/vRupture;

tslip = 0.0*Units.hr; % time of slip starting
t0 = 0.0; % time to record temperature

% interval of temperature recording
dt = [1*Units.s, 1*Units.min, 1*Units.hr, 1*Units.day]; 

% number of temp recordings
nt = [60, 59*48, 24*7, 365]; 

%% program
fprintf( 'Program %s\n', mfilename );

%Calculate the heat generated by slip
heat = stress*slip*dx*dz;
fprintf( 'Heat from slip = %.4g\n', heat);

%Calculate the heat rate for the eqk
tEndEqk = tslip + eqkSlipTime;
dqdt_eqk = heat/eqkSlipTime;
fprintf('Heat rate from eqk = %.4g\n', dqdt_eqk);

%Calculate heat rate for the creep
tEndCreep = tslip + creepSlipTime;
dqdt_creep = heat/creepSlipTime;
fprintf('Heat rate from creep = %.4g\n', dqdt_creep );

%Calculate heat rate for the second creep
tEndCreep2 = tslip + creepSlipTime2;
dqdt_creep2 = heat/creepSlipTime2;
fprintf('Heat rate from creep 2 = %.4g\n', dqdt_creep2 );

% Loop through all measurement points
nOut = sum( nt );
tout = zeros(nOut,1);
out = zeros(nOut, 3);
iOut = 0;

for j=1:4,
    for i=0:nt(j)-1, 
        % Current time
        iOut = iOut+1;
        t = t0 + i*dt(j);
            
        % calc temperature
        dT_eqk  = tempfromheat( dqdt_eqk, dx*dz, tslip, tEndEqk, t, ...
            faultWidth, diffusivity, density, specHeat );
        dT_creep  = tempfromheat( dqdt_creep, dx*dz, tslip, tEndCreep, t,  ...
            faultWidth, diffusivity, density, specHeat ); 
        dT_creep2  = tempfromheat( dqdt_creep2, dx*dz, tslip, tEndCreep2, t,  ...
            faultWidth, diffusivity, density, specHeat ); 
        
        %fprintf('%8.4f, %.3f, %.3f, %3f\n', t/Units.hr, dT_eqk, dT_creep, dT_creep2);
        tout(iOut) = t;
        out(iOut,:) = [dT_eqk, dT_creep, dT_creep2];
    end
    t0 = t;
end
  
%%
figure
plot( tout./Units.s, out, '-+','linewidth', 2 )
set( gca, 'xscale', 'log')
xlabel('Time [s]')
ylabel('Temperature Change [k]')
grid on;
return
% test_cooling.m ends here