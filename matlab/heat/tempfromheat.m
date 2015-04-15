function temp = tempfromheat( heatRate, cellArea, t1, t2, tNow, faultWidth, ...
    diffusivity, density, specHeat )
%
% temp = tempfromheat( heatRate, cellArea, t1, t2, tNow, faultWidth, diffusivity )
%

temp = inst_temp(heatRate, cellArea, faultWidth, density, specHeat)*...
    integ_tempdecay( t1, t2, tNow, faultWidth, diffusivity );
  
end
% ------------------------------------------------------------------------------
function tempSpike = inst_temp(heatRate, cellArea, faultWidth, density, specHeat)
% Instantaneous temperature from heat rate
tempSpike = 0.5*heatRate./( cellArea*faultWidth*density*specHeat );

end

% ------------------------------------------------------------------------------
function result = integ_tempdecay( tSlipStart, tSlipEnd, tNow, faultWidth, diffusivity )
% Check what time to return heat for */

% Error tolerance
EPS = 1e-12;
 
if( tNow <= tSlipStart),
    % Before t=0
    result = 0.0;
else
    % After t=0
    tSinceStart = tNow - tSlipStart;
    tSinceEnd = max(0.0, tNow - tSlipEnd);

    if( abs( tSinceEnd-tSinceStart ) < EPS*tSinceEnd ),
      % Slip pulse
      result = (tSinceStart-tSinceEnd)*tempDecay( 0.5*(tSinceStart+tSinceEnd) );
    else
      % Period of slip
      result = integral( @(x) tempdecay(x, faultWidth, diffusivity), ...
          tSinceEnd, tSinceStart, 'RelTol', EPS );
    end
end

end
%------------------------------------------------------------------------------
function result = tempdecay( tSinceSlip, faultWidth, diffusivity )
% Temperature decay function
  
% Check greater than 0.0 
if( tSinceSlip > 0.0 ),
    tmp = 0.25*faultWidth./sqrt( diffusivity*tSinceSlip );
    result = erf(tmp) - erf(-tmp);
else
    % Case where tSinceSlip = 0.0: erf(Inf) - erf(Inf) = 1 - -1 = 2
    result = 2.0;
end

end
