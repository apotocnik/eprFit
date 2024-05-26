%--------------------------------------------------------------------------
% markers  - return markers simbol acording to number i
%
% Author: Anton Potocnik, F5, IJS
% Date:   15.03.2009
% Arguments:
%       m = markers(i)
% Input:    i       ... integer 1 ... Inf
% Output:   m       ... marker's simbol (+,o,*,x,...)
%--------------------------------------------------------------------------


function m = markers(i)

i = mod(i,13)+1;
mar = {'o','^','s','v','d','>','<','+','*','x','p','h','.'};
m = mar{i};


