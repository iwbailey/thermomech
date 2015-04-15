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
    's', 1, ...
    'GPa', 1e9, ...
    'MPa', 1e6, ...
    'km', 1e3, ...
    'mm', 1e-3, ...
    'kJ', 1e3 );

%% Parameters used in this program
nx = 128;
nz = 32;
faultLength = 70.0*Units.km;
faultHeight = 17.5*Units.km;
faultWidth = 100*Units.mm;
surfTemp = 273 + 20;
dTdz = 20 *(1/Units.km);
rigidity = 30*Units.GPa;
tau0 = 0.6*Units.MPa; % cohesion i.e. static shear strength at surface
dSigmadz = 18*(Units.MPa/Units.km); % Pa/m
ifile_strengthdrops = fullfile('..','inputs','stressdrops_unif.txt');
zBD = 7.5*Units.km; % Brittle-ductile length [m], TODO: check
fs = 0.75; % coefficient of friction
dosCoef = 1.25; % Dynamic overshoot coefficient
faultE = 0.0; % No thermal impact, set activation energy to zero
faultn = 3;
plateVelocity = 35 *( Units.mm/Units.yr );
slipVelocity = 6 *(Units.km/Units.s);
nTimeMax = 1e6;
time0 = 125 *Units.yr;
maxTime = 50 *Units.yr;
minTimeStep = 1 *Units.min;
maxTimeStep = 1 *Units.yr;

%% Initialization
% Print program name
fprintf('%s:\n', mfilenamename)

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
strengthDrops = reshape( tmp(:,3), nx, nz)';
strengthDynamic = strengthStatic - strengthDrops;

% Compute the stiffness matrix
[stiffk, selfstiff] = stiffnessmatrixk( cellLength, cellHeight, depths(1), ...
    rigidity, nx, nz);
% [stiff, selfstiff1] = stiffnessmatrix( cellLength, cellHeight, depths(1), ...
%     rigidity, nx, nz);

% Calculate the A value according to the  creep mask used in BZ1996
faultA = faultcreep_bz96( nx, nz, cellLength, cellHeight, ...
                          zBD, zBD, plateVelocity, ...
                          tau0 + fs*dSigmadz*(faultHeight-zBD) );

% Calculate the creep strength at the loading strain rate and bkgd temp
strengthCreep = creepstrength( plateVelocity, faultA, faultn, faultE, bkgdT );

% Set the initial stress
initStress = min( 0.9*strengthStatic, 0.95*strengthCreep );

% Set containers for recording
totalCreepSlip = zeros( nz, nx );
totalEqkSlip = zeros( nz, nx);

%% Initialize the fault
% Set up structure for the static properties of the fault
Flt = struct( ...
    'nL', nx, 'nD', nz, ...
    'cellLength', cellLength, 'cellHeight', cellHeight, 'width', faultWidth, ...
    'strengthStatic', strengthStatic, 'strengthDynamic', strengthDynamic, ...
    'stiffness', selfstiff, 'dynOvershootCoeff', dosCoef, ...
    'arrhA', faultA, 'stressExpon', faultn, 'activEnergy', faultE, ...
    'initStress', initStress, 'initTemp', bkgdT );

% Dynamic properties of the fault
Flt.creepVel = creeprate( Flt.initStress, Flt.arrhA, Flt.stressExpon, ...
    Flt.activEnergy, Flt.initTemp );
Flt.slipDeficit = time0*(plateVelocity*ones(nz,nx)-Flt.creepVel);
Flt.stress = Flt.initStress + slipdeftostressk( Flt.slipDeficit, stiffk, nx, nz );
%Flt.stress2 = Flt.initStress + slipdeftostress( Flt.slipDeficit, stiff, nx, nz );

%% Start the algorithm

% Print header for the earthquake catalog
fprintf('Time_yr, x_km, z_km, Mag_P, Mag_W, Area_km2, StressDrop_MPa\n');

% Run the algorithm
iTimeStep = 0;
time = 0.0;
while( iTimeStep < nTimeMax && time < maxTime )

    % Calculate whether there are any hypocenters
    nCrit = nnz( Flt.stress >= Flt.strengthStatic );
    if( nCrit > 1 )
        fprintf('\nWARNING: Multiple hypocenters (%i)\n', nCrit );
    end
    
    if( nCrit > 0 ),
        % Output  progress to terminal
        %fprintf('*EQ* ');

        % Hypocenter
        overStress = Flt.stress - Flt.strengthStatic;
        [iHypo, jHypo] = find( overStress == max( overStress(:) ) );
        
        % Compute the slip deficit before the earthquake
        slipDeficitPre = Flt.slipDeficit;
        stressPre = Flt.stress;
        timePre = time;

        % Compute the earthquake
        %eqkSlip = calcearthquake( Flt, stiffk, time, slipVelocity, stiff );
        eqkSlip = calcearthquake( Flt, stiffk, time, slipVelocity );


        % Compute the new slip deficit
        slipDeficitPost = slipDeficitPre - eqkSlip;
        
        Flt.slipDeficit = slipDeficitPost;
%         Flt.stress =  Flt.initStress + slipdeftostressk( Flt.slipDeficit, ...
%                                                        stiffk, nx, nz );
        Flt.stress =  Flt.initStress + slipdeftostressk( Flt.slipDeficit, ...
            stiffk, nx, nz );
        
        %         % Record the total earthquake slip
        %         if( iTimeStep > 0 )
        %             totalEqkSlip = totalEqkSlip + (slipDeficitPre - slipDeficitPost );

        % Output the earthquake catalog: time, x, z, magP, magW, area, stressdrop*/
        fprintf( '%9.6f, %6.2f, %5.2f, %4.2f, %4.2f, %6.2f, %5.2f\n',...
            timePre/Units.yr,...
            (jHypo-0.5)*Flt.cellLength/Units.km ,...
            (iHypo-0.5)*Flt.cellHeight/Units.km ,...
            potmagnitude( scalarpotency( slipDeficitPre, slipDeficitPost, Flt.cellLength*Flt.cellHeight) ), ...
            mommagnitude( scalarmoment( slipDeficitPre, slipDeficitPost, Flt.cellLength*Flt.cellHeight, rigidity) ), ...
            Flt.cellLength*Flt.cellHeight*nnz(eqkSlip>0) / (Units.km*Units.km), ...
            mean( stressPre(eqkSlip>0) - Flt.stress(eqkSlip>0) )/Units.MPa );
    end %END earthquake
    
    % Get the creep velocity on the fault
    Flt.creepVel = Flt.width*creeprate( Flt.stress, Flt.arrhA, Flt.stressExpon, ...
                                          Flt.activEnergy, Flt.initTemp );

    % Compute the time to failure
    timeStep = timetofailure( Flt, Flt.stress, Flt.creepVel, plateVelocity, ...
                              stiffk );

    % Adjust the time step
    if( timeStep < 0 )
        % Negative implies something went wrong
        error('Negative time step at iTimeStep %i, t = %.3f yr\n', ...
              iTimeStep, time / Units.yr );
    elseif( timeStep > maxTimeStep )
        % Don't let it get too big or creep rates will be inaccurate
        %fprintf('Time %8.2f: Using maximum time step\n', time/Units.yr );
        timeStep = maxTimeStep;
    elseif( timeStep < minTimeStep )
        % Don't let it get too small or we will be waiting
        %fprintf('Time %8.2f: Using minimum time step\n', time/Units.yr );
        timeStep = minTimeStep;
    %else
        %fprintf('Time %8.2f\n', time/Units.yr);
    end

    % Load the fault
    Flt.slipDeficit = Flt.slipDeficit + ...
                         timeStep*(plateVelocity - Flt.creepVel);

    % Add on the creep slip
    totalCreepSlip = totalCreepSlip + Flt.creepVel*timeStep;

    % Update the time
    time = time+timeStep;

    % Get the new stress
    Flt.stress =  Flt.initStress + slipdeftostressk( Flt.slipDeficit, ...
        stiffk, nx, nz );

    iTimeStep = iTimeStep + 1;

%     %fprintf('\n');
%     % TMP
%     clf
%     subplot(2,1,1);
%     imagesc( Flt.stress-Flt.initStress ); axis equal tight;
%     cb = colorbar;
%     title('Stress - Init Stress')
%     ylabel(cb, 'MPa')
%     
%     subplot(2,1,2);
%     imagesc( Flt.slipDeficit ); axis equal tight;
%     cb = colorbar;
%     title('SlipDeficit');
%     ylabel(cb, '[mm]')
%     pause(0.1)

end

% Break the line
fprintf('\n');

% Report why we finished the loop
if( iTimeStep == nTimeMax ), fprintf('Max number of iterations reached\n');
else fprintf( 'Max time reached\n');
end

%% Write Output files

% Total creep slip
ofile_creepslip = sprintf('%s_%s.txt', mfilename, 'creepslip' );
fid = fopen( ofile_creepslip, 'w');
fprintf( fid, '%.5e\n', totalCreepSlip(:)' );
fclose(fid);

% Total earthquake slip
ofile_eqslip =  sprintf('%s_%s.txt', mfilename, 'eqslip' );
fid = fopen( ofile_eqslip, 'w');
fprintf( fid, '%.5e\n', totalEqkSlip(:));
fclose( fid );

% Remaining slip deficit
ofile_slipdef =  sprintf('%s_%s.txt', mfilename, 'slipdef' );
fid = fopen( ofile_slipdef, 'w');
fprintf( fid, '%.5e\n', Flt.slipDeficit(:) );
fclose( fid );

% Current stress
ofile_stress =  sprintf('%s_%s.txt', mfilename, 'stress_MPa' );
fid = fopen( ofile_stress, 'w');
fprintf( fid, '%.5e\n', Flt.stress(:)./Units.MPa );
fclose( fid );

% END
