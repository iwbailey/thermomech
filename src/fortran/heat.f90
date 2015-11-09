module heat
  ! Error tolerance
  real(kind=8), parameter :: EPS = 1e-12;

contains
  !-------------------------------------------------------------------------------
  function integ_heatdecay(SlipStart, tSlipEnd, tNow, faultWidth, diffusivity )
    !
    !
    ! result = integ_tempdecay( tSlipStart, tSlipEnd, tNow, faultWidth, diffusivity )

    ! Check what time to return heat for

    ! Check tStart and tEnd for errors
    if( tSlipEnd<=tSlipStart) then
       tSlipEnd = tSlipStart + EPS;
       warning('Changed tSlipEnd to %f', tSlipEnd)
    end if

    if( tNow <= tSlipStart ) then
       integ_heatdecay = 0;
    else
       ! Get the time since slip
       tSinceStart = tNow - tSlipStart;

       ! Get the time since the end of slip. If not yet finished set to zero
       tSinceEnd = max(0.0, tNow - tSlipEnd);

       ! Find cases where the slip time is much smaller than the time since finished
       slipTime = tSinceStart - tSinceEnd;

       if( slipTime < EPS*tSinceEnd ) then
          ! slip pulse
          integ_heatdecay = slipTime* &
               tempdecay( tSinceEnd, faultWidth, diffusivity );
       else
          ! Period of slip
          integ_heatdecay = quadcc( @(x) tempdecay(x, faultWidth, diffusivity), &
               tSinceEnd, tSinceStart, EPS );
       end if

    end if
  end function integ_heatdecay

  !------------------------------------------------------------------------------
  function tempdecay( tSinceSlip, faultWidth, diffusivity )
    ! Temperature decay function

    ! Case where tSinceSlip = 0.0: erf(Inf) - erf(Inf) = 1 - -1 = 2
    result = 2*ones(size(tSinceSlip));

    ! Check greater than 0.0
    isOk = tSinceSlip > 0;
    tmp = 0.25*faultWidth./sqrt( diffusivity*tSinceSlip(isOk) );
    result(isOk) = erf(tmp) - erf(-tmp);

  end function tempdecay


end module heat
