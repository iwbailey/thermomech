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
