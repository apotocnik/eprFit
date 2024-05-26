%--------------------------------------------------------------------------
% xc2g  - calculates g from frequency and xc (Hc) (H)
%
% Author: Anton Potocnik, F5, IJS
% Date:   27.01.2009 - 22.03.2010
% Arguments:
%       g = xc2g(xc,freq)
%           xc in mT!!!
%           freq in GHz
%       if freq == 0, g = xc;
%--------------------------------------------------------------------------

function g = xc2g(xc,freq)

    if freq > 0
        g = 71.44775*freq./xc;
    else
        g = xc;
    end
