function result = integ_tempdecay( tSlipStart, tSlipEnd, tNow, faultWidth, diffusivity )
%INTEG_TEMPDECAY
%
% result = integ_tempdecay( tSlipStart, tSlipEnd, tNow, faultWidth, diffusivity )

% Check what time to return heat for

% Error tolerance
EPS = 1e-12;

result = zeros(size(tNow));
isNZ = tNow>tSlipStart;

if( ~any(isNZ) ), return; end

% Check tStart and tEnd for errors
if( tSlipEnd<=tSlipStart),
    tSlipEnd = tSlipStart + EPS;
    warning('Changed tSlipEnd to %f', tSlipEnd)
end

% Get the time since slip
tSinceStart = tNow - tSlipStart;

% Get the time since the end of slip. If not yet finished set to zero
tSinceEnd = max(0.0, tNow - tSlipEnd);

% Find cases where the slip time is much smaller than the time since finished
slipTime = tSinceStart - tSinceEnd;
isSlipPulse = slipTime < EPS*tSinceEnd;

% Get those results
result(isSlipPulse) = slipTime(isSlipPulse).*...
    tempdecay( tSinceEnd(isSlipPulse), faultWidth, diffusivity );

idx = find( isNZ & ~isSlipPulse)';
for i1=1:numel(idx),
    i = idx(i1);
    % Period of slip
    %      result = integral( @(x) tempdecay(x, faultWidth, diffusivity), ...
    %    tSinceEnd, tSinceStart, 'RelTol', EPS );
    result(i) = quadcc( @(x) tempdecay(x, faultWidth, diffusivity), ...
                        tSinceEnd(i), tSinceStart(i), EPS );

end

end
%------------------------------------------------------------------------------
function result = tempdecay( tSinceSlip, faultWidth, diffusivity )
% Temperature decay function

% Case where tSinceSlip = 0.0: erf(Inf) - erf(Inf) = 1 - -1 = 2
result = 2*ones(size(tSinceSlip));

% Check greater than 0.0
isOk = tSinceSlip > 0;
tmp = 0.25*faultWidth./sqrt( diffusivity*tSinceSlip(isOk) );
result(isOk) = erf(tmp) - erf(-tmp);

end
