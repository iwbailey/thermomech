%
%---  test_bz1996.cpp ---
%
%   Test a fault set up that is analogous to the BZ 1996 model,
%   i.e. stress-dependent, but no temperature dependent creep.  Note it is not
%   exactly the same because the creep is set up in a different way
%
%
%  Filename: script_bz1996.m
%  Description:
%  Author: Iain W. Bailey
%  Created: Sun Feb  18 20:39:13 2015 (-0500)
%  Version: 1
%  Last-Updated: Sun Mar  8 17:46:42 2015 (-0400)
%            By: Iain W. Bailey
%      Update #: 297
%

%  Change Log:
%
%
%
%

%% Define units used here
Units = struct(...
    'yr', 365.25*24*3600,...
    'day', 24*3600, ...
    'min', 60, ...
    'GPa', 1e9, ...
    'MPa', 1e6, ...
    'km', 1e3, ...
    'mm', 1e-3, ...
    'kJ', 1e3 );

%% Parameters used in this program
surfTemp = 273 + 20;
dTdz = 20 *1/Units.km;
nx = 128;
nz = 32;
faultLength = 70.0*Units.km;
faultHeight = 17.5*Units.km;
faultWidth = 100*Units.mm;
rigidity = 30*Units.GPa;
tau0 = 6*Units.MPa; % cohesion i.e. static shear strength at surface
dSigmadz = 18*(Units.MPa/Units.km); % Pa/m
ifile_strengthdrops = fullfile('..','inputs','stressdrops_unif.txt');
zBD = 7.5*Units.km; % Brittle-ductile length [m], TODO: check
fs = 0.75; % coefficient of friction
dosCoef = 1.25; % Dynamic overshoot coefficient
faultE = 0.0; % No thermal impact, set activation energy to zero
faultn = 3;
plateVelocity = 35*( Units.mm/Units.yr );
nTimeMax = 1e6;
maxTime = 1e5 * Units.yr;
minTimeStep = 1 *Units.min;
maxTimeStep = 1 *Units.day;

%% Initialization
% Print program name
fprintf('%s:\n', mfilename)

% Depth of mid point of each cell
cellHeight = faultHeight/nz;
cellLength = faultLength/nx;
depths = (faultHeight/nz)* (0.5:1:nz-0.5)';

% Background temperature
bkgdT = repmat( surfTemp + dTdz*depths, 1, nx);

% Static strength profile
strengthStatic = repmat( tau0 + fs*dSigmadz*depths, 1, nx);

% Dynamic strength drop
tmp = load( ifile_strengthdrops );
strengthDrops = resize(tmp(:,3), nx, nz)';
strengthDynamic = strengthStatic - strengthDrops;

% Compute the stiffness matrix
stiffnessMatrix = stiffnessmatrix( cellLength, cellHeight, depths(1), ...
                                  rigidity, nx, nz);


% Calculate the A value according to the  creep mask used in BZ1996
faultA = faultcreep_bz96( nx, nz, cellLength, cellHeight, ...
                          zBD, zBD, plateVelocity, ...
                          tau0 + fs*dSigmadz*(faultHeight-zBD) );

% Set up the creep parameters over entire fault including creep mask

% Calculate the creep strength at the loading strain rate and bkgd temp
strengthCreep = creepstrength( plateVelocity, faultA, faultn, faultE, bkgdT );

% Set the initial stress
initStress = min( strengthStatic, strengthCreep );

%% Initialize the fault
% Set up structure for the static properties of the fault
F0 = struct( ...
    'nL', nx, 'nD', nz, ...
    'cellLength', cellLength, 'cellHeight', cellHeight, 'width', faultWidth, ...
    'strengthStatic', strengthStatic, 'strengthDynamic', strengthDynamic, ...
    'dynOvershootCoeff', dosCoef, ...
    'arrhA', faultA, 'stressExpon', faultn, 'activEnergy', faultE, ...
    'initStress', initStress, 'initTemp', bkgdT );

% Dynamic properties of the fault
FltDyn = struct( ...
    'slipDeficit',zeros(nz, nx), ...
    'stress', F0.initStress + ...
    slipdeftostress( zeros(nz,nx), stiffnessMatrix, nx, nz ) );

%% Set containers for recording
totalCreepSlip = zeros( nz, nx );
totalEqkSlip = zeros( nz, nx);

%% Start the algorithm

% Print header for the earthquake catalog
fprintf('Time_yr, x_km, z_km, Mag_P, Mag_W, Area_km2, StressDrop_MPa\n');

% Run the algorithm
iTimeStep = 0;
time = 0.0;
while( iTimeStep < nTimeMax && time < maxTime )

    % Get the creep velocity on the fault
    FltDyn.creepVel = F0.width*creeprate( FltDyn.stress, F0.arrhA, F0.stressExpon, ...
                                          F0.activEnergy, F0.initTemp );

    % Compute the time to failure
    timeStep = timetofailure( F0, FltDyn.stress, FltDyn.creepVel, plateVelocity, ...
                              stiffnessMatrix );

    % Adjust the time step
    if( timeStep < 0 )
        % Negative implies something went wrong
        error('Negative time step at iTimeStep %i, t = %.3f yr\n', ...
              iTimeStep, time / Units.yr );
    elseif( timeStep > maxTimeStep )
        % Don't let it get too big or creep rates will be inaccurate
        fprintf('Time %8.2f: Using maximum time step\n', time/Units.yr );
        timeStep = maxTimeStep;
    elseif( timeStep < minTimeStep )
        % Don't let it get too small or we will be waiting
        fprintf('Time %8.2f: Using minimum time step\n', time/Units.yr );
        timeStep = minTimeStep;
    else
        fprintf('Time %8.2f\n', time/Units.yr);
    end

    % Load the fault
    FltDyn.slipDeficit = FltDyn.slipDeficit + (plateVelocity*timeStep - ...
                                               FltDyn.creepVel);

    % Add on the creep slip
    totalCreepSlip = totalCreepSlip + FltDyn.creepVel*timeStep;

    % Update the time
    time = time+timeStep;

    % Get the new stress
    FltDyn.stress =  F0.initStress + slipdeftostress( FltDyn.slipDeficit, ...
                                                      stiffnessMatrix, nx, nz ...
                                                      );

    % Calculate whether there are any hypocenters
    nCrit = nnz( FltDyn.stress >= F0.strengthStatic );
    if( nCrit > 1 )
        fprintf('\nWARNING: Multiple hypocenters (%i)\m', nCrit );
    end

    if( nCrit > 0 )
        % Output  progress to terminal */
        fprintf('*EQ* ');

        % Compute the slip deficit before the earthquake */
        slipDeficitPre = FltDyn.slipDeficit;
        stressPre = FltDyn.stress;
        timePre = time;

        % Compute the earthquake */
        eqkSlip = computeEarthquake( stress, time, iHypo, jHypo, slipVelocity );

        % Compute the new slip deficit */
        slipDeficitPost = slipDeficitPre - eqkSlip;

        % Compute the earthquake properties */
        thisEQ = Earthquake( iHypo, jHypo, timePre, faultLength*faultWidth/nx/nz, ...
                             slipDeficitPre, slipDeficitPost, stressPre, stress );

        % Record the total earthquake slip */
        if( iTimeStep > 0 )
            totalEqkSlip = totalEqkSlip + (slipDeficitPre - slipDeficitPost );

            % Output the earthquake catalog: time, x, z, magP, magW, area, stressdrop*/
            fprintf( '%9.6E, %6.2f, %5.2f, %4.2f, %4.2f, %6.2f, %5.2f\n',
                  timePre/Units.yr,
                  (iHypo+0.5)*cellLength/Units.km ,
                  (jHypo+0.5)*cellHeight/Units.km ,
                  thisEQ.potMagnitude() ,
                  thisEQ.momMagnitude( rigidity ) ,
                  thisEQ.ruptureArea() / (Units.km*Units.km) ,
                  thisEQ.staticStressDrop()/Units.MPa );
        end
      end %END earthquake

      iTimeStep = iTimeStep + 1;

      fprintf('\n');

end

% Break the line
fprintf('\n');

% Report why we finished the loop */
if( i == nTimeMax ), fprintf('Max number of iterations reached\n');
else fprintf( 'Max time reached\n');
end

%% Write Output files

% Total creep slip
ofile_creepslip = sprintf('%s_%s.txt', mfile, ofilesuffix_creepSlip );
fid = fopen( ofile_creepslip, 'w');
fprintf( fid, '%.5e\n', totalCreepSlip(:)' );
fclose(fid);

% Total earthquake slip
ofile_eqslip =  sprintf('%s_%s.txt', mfile, ofilesuffix_eqSlip );
fid = fopen( ofile_eqslip, 'w');
fprintf( fid, '%.5e\n', totalEqkSlip(:));
fclose( fid );

% Remaining slip deficit
ofile_slipdef =  sprintf('%s_%s.txt', mfile, ofilesuffix_slipdef );
fid = fopen( ofile_slipdef, 'w');
fprintf( fid, '%.5e\n', Fault.slipDeficit(:)./Units.m );
fclose( fid );

% Current stress
ofile_stress =  sprintf('%s_%s.txt', mfile, ofilesuffix_stress );
fid = fopen( ofile_stress, 'w');
fprintf( fid, '%.5e\n', Fault.stress(:)./Units.MPa );
fclose( fid );

% END
