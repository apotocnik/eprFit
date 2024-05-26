%--------------------------------------------------------------------------
% units2value  - Converts string units to its numerical values
%
% Version: 1.0
% Author: Anton Potocnik, F5, IJS
% Date:   25.01.2009
% Arguments:
%       units2value(txt)
%--------------------------------------------------------------------------
function Y = units2value(txt)

switch txt
    case 'n'
        Y = 1e-9;
    case 'u'
        Y = 1e-6;
    case 'm'
        Y = 1e-3;
    case 'k'
        Y = 1e3;
    case 'M'
        Y = 1e6;
    case 'G'
        Y = 1e9;
    otherwise
        Y = 1;
end