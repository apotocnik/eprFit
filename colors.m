%--------------------------------------------------------------------------
% colors  - return colors simbol acording to number i
%
% Author: Anton Potocnik, F5, IJS
% Date:   15.03.2009
% Arguments:
%       c = colors(i)
% Input:    i       ... integer 1 ... Inf
% Output:   c       ... colors's simbol (b,g,r,c,...)
%--------------------------------------------------------------------------


function c = colors(i)

colo = 'bkgrcmy';
c = colo(mod(i,7)+1);


