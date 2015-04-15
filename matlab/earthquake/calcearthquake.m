function [totalSlip, time] = calcearthquake( Flt, stiffk, time, slipVelocity )
% Compute the cascading failure of cells given an initial "hypocenter" cell
% already in a critical state.  Update the slip deficit, stress and time.
% Output the distribution of total slip on the fault.
%
% IN/OUT
%  stress : current stress on fault, updated by function
%  time : current time, updated according to the max amount of slip and the slip
%  velocity
% IN
%  iHypo : along strike index of the hypocenter cell
%  jHypo : down-dip index of the hypocenter cell
%  slipVelocity: constant slip velocity of seismic slip in the model
% OUT
%  totalSlip : slip in each cell of the fault during the earthquake
%
maxNit = 1e5;

% Check that the hypocenter is critical
isCritical = ( Flt.stress >= Flt.strengthStatic );

% Container for recording the total slip
totalSlip = zeros( Flt.nD, Flt.nL );

% Set up the first critical cell 
isStatic = true( Flt.nD, Flt.nL );

arrestStress = Flt.strengthStatic - ...
    Flt.dynOvershootCoeff*(Flt.strengthStatic - Flt.strengthDynamic);

% Compute stress redistribution and slip until all cells are below
% strength
for i=1:maxNit, 
    if( ~any(isCritical) ), break; end
    
%     % TMP
%     subplot(3,1,1);
%     imagesc( Flt.stress ); axis equal tight;
%     cb = colorbar;
%     subplot(3,1,2);
%     imagesc( Flt.slipDeficit ); axis equal tight;
%     cb = colorbar;
%     subplot(3,1,3);
%     imagesc( Flt.stress - Flt.stress2 ); axis equal tight;
%     cb = colorbar;

    % Get the stress drop of failed cells
    stressDrops = Flt.stress(isCritical) - arrestStress(isCritical);
    
    % Compute the amount of slip due to the stress drop
    slip = stressDrops ./ -Flt.stiffness(isCritical);
    
    % Failed cells now have dynamic strength
    isStatic(isCritical) = false;
    
    % Adjust the slip deficit and total slip
    Flt.slipDeficit(isCritical) = Flt.slipDeficit(isCritical) - slip;
    totalSlip(isCritical) = totalSlip(isCritical) + slip;

    % Recompute the stress 
    Flt.stress = Flt.initStress + slipdeftostressk( Flt.slipDeficit, stiffk );
    %Flt.stress2 = Flt.initStress + slipdeftostress( Flt.slipDeficit, stiff );

    if( any(Flt.stress <0 )), error('Negative Stress'); end
    
    % Find any failed cells 
    isCritical = (isStatic & Flt.stress >= Flt.strengthStatic) | ...
        (~isStatic & Flt.stress >= Flt.strengthDynamic);
end
if( i == maxNit ), 
    error('Maximum number of iterations reached');
end
%     % TMP
%     subplot(2,1,1);
%     imagesc( Flt.stress ); axis equal tight;
%     cb = colorbar;
%     subplot(2,1,2);
%     imagesc( Flt.slipDeficit ); axis equal tight;
%     cb = colorbar;

% Update time based on the maximum slip amount
time = time + max( totalSlip(:) )/slipVelocity;

end