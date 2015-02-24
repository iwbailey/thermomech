function time = timetofailure( Fault, stress, creepVel, plateVelocity, stiffnessMatrix )
% Estimate time to failure for the entire fault given loading at current rate
% and stress. Note that this does not account for changes to the creep velocity
% caused by changes to the stress. Hence, this function should not be used to
% set the time step unless constrained by an upper bound.

% Get amount of stress before failure
stressToFailure = Fault.strengthStatic - stress;

% Check if already at failure
if( min( stressToFailure(:) ) <= 0 ), time = 0.0; return; end

% Convert the loading rate from all cells on the fault into a stress build up
% rate for this fault, assuming constant loading
loadingRate = plateVelocity - creepVel;

stressingRate = slipdeftostress( -loadingRate, stiffnessMatrix, Fault.nL, Fault.nD );

% Compute the time taken to reach failure
failTime = stressToFailure ./ stressingRate;

% Get min positive time
time = min( failTime( stressingRate>0 ));

end
