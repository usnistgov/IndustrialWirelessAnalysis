%
%   Compute the RFNest D508 profile file
%   distance offset in meters must be under 1 km
%
%
% Author: Rick Candell
% National Institute of Standards and Technology

% tap delay Doppler spread distribution

NTAPS = 13;
CHAN_TYPE = 2;  % 0 for freespace, 1 rayleigh, 2 rician

delay_profile = ...
    [   
        CHAN_TYPE zeros(1,NTAPS-1);
        floor(32767 * (10.^([-12 0 -27 -6 -15 -10 -18 -10 -24 -18 -35 -12 -30]./20)));
        round(65535 * rand(1, NTAPS));
        zeros(1,NTAPS);
        [40, 160, 200, 225, 300, 310, 350, 380, 400, 425, 480, 580, 800] ;
    ]

csvwrite('rician_cir01.csv', delay_profile(:,:))
