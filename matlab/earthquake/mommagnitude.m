function mag = mommagnitude( moment )
%MOMMAGNITUDE Moment magnitude, given the moment
% log10 M0 = 1.5 mW + 16.1 
%
% mag = mommagnitude( moment )

logM = log10( moment );
mag = (logM - 9.1)/1.5;
end
