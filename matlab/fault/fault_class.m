classdef fault_class
%
% Class definition for Fault object.  Holds grid of slip cells, their
%    static and evolving properties.  Base class of ThermoFault
%   A fault is a grid of slip surfaces with no dip embedded in an
%   elastic halfspace
%
properties ( SetAccess=protected,GetAccess=public)
    nL = 128;
    nD = 32;
    dx = 70e3/nL;
    dz = 17.5e3/nD;
    faultWidth = 10e-2;
    dynamicOvershootCoeff = 1.25;
    initStress = zeros(nD, nL); % initial stress for each cell
    strengthStatic = zeros( nD, nL); % static strength
    strengthDynamic = zeros( nD, nL); % dynamic strength
    arrhA = zeros( nD, nL); % Arrhenius amplitude
    stressExpon = 3*ones(nD,nL); % exponent to stress term in arrhenius relation
    activEnergy = 130e3*ones(nD,nL); % Activation energy
    slipDeficit = zeros(nD,nL);
    stiffnessMatrix;

end

properties (Dependent)
    stress
    creepVelocity
    arrestStress
end

methods
    %----------------------------------------------------------------------
        function obj = fault_class( varargin )
            % Constructor
            if(nargin==1 && isstruct(varargin{1}))
                % Initialised from a structure
                for f = fieldnames( varargin{1} )',
                    fn = char(f);
                    obj.(fn) = varargin{1}.(fn);
                end
            end
            % ...
        %     nL, ...% along-strike number of grid cells
        %     nD, ...% down-dip number of grid cells
        %     dx, ...% length of a single slip cell
        %     dz, ...% height of a single slip cell
        %     w, ...% fault width for strain rate calc
        %     staticStrengths, ...% Static strength for each slip cell
        %     dynamicStrengths, ...% dynamic strength for each slip cell
        %     dynamicOvershootCoeff, ...% dynamic overshoot coefficient for fault
        %     arrheniusAmpl, ...% amplitude in the arrhenius relation
        %     stressExpon, ...% Stress exponent in the arrhenius relation
        %     actEnergy, ...% activation energy
        %     stiffnessMatrix, ...% Stiffness matrix for the fault
        %     initStress) % initial temperature for each slip cell
        end

    %----------------------------------------------------------------------
        function loadFault( obj, plateVelocity, creepVel, dtime )
        % Load the fault with constant plate motion for a fixed time period when creep
        %    velocity has already been calculated
        obj.slipDeficit = obj.slipDeficit + dtime*(plateVelocity-creepVel);
        end

    %----------------------------------------------------------------------
        function stress = get.stress( obj )
        % Get the stress for each slip cell on the fault

        % TODO
        stress = zeros(obj.nD, obj.nL);
        end

    %----------------------------------------------------------------------
        function creepVel = get.creepVelocity( obj, stress, temperature )
        % Get the current creep velocity on the fault

        % TODO
        stress = zeros(obj.nD, obj.nL);
        end
    %----------------------------------------------------------------------
        function taua = get.arrestStress( obj )
        % Get the arrest stress
        taua = (obj.strengthStatic - obj.strengthDynamic)*obj.dynamicOvershootCoeff;
        end

    %----------------------------------------------------------------------
        function time =  estimateTimeToFailure( obj, stress, creepVel, plateVelocity)
        % Estimate time to failure for entire fault

        % TODO
        time = 1;
        end

        function idx = findCriticalCells( obj, stress, isStatic )
        % Find cells with stress above their failure strength

        % Get the current stress
        faultStress = obj.stress;

        % Find the cells where stress is above the appropriate failure strength
        idx = find( (isStatic & faultStress >= strengthStatic) | ...
                    (~isStatic & faultStress >= strengthDynamic) );
        end

        function [totalSlip, time, stress] = computeEarthquake( obj, ...
                                                          stress, time, iHypo, ...
                                                          jHypo, ...
                                                          slipVelocity)
        % Compute an earthquake on the fault once hypocenter(s) are known

        % TODO
        end
end


end
%------------------------------------------------------------------------------
