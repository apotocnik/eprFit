%--------------------------------------------------------------------------
% extrange  -  Extract range from data
%
% Author: Anton Potocnik, F5, IJS
% Date:   01.10.2009
% Arguments:
%       [H,Y,mask] = extrange(H, Y, range])
% Input:
%        H, Y   - specter
%        range  - string of ranges "[100 120] [500 550] [1000 3000]"
% Output:
%        H, Y  - truncated specter
%        mask  - [0,0,1,1,1,1,0,1,1,1,1] mask which are in range
%--------------------------------------------------------------------------

function [H Y mask] = extrange(H,Y,range)

    range = eval(['[' range ']']);
    mask = [];
	for i=1:2:numel(range)
        if range(i) == -Inf
            [C,idx_fst] = min(H);  % get lowest first
        else
            [C,idx_fst] = min(abs(H-range(i)));  % get nearest first
        end
        if range(i+1) == Inf
            [C,idx_lst] = max(H);  % get higest last
        else
            [C,idx_lst] = min(abs(H-range(i+1)));  % get nearest last
        end
        mask = [mask idx_fst:idx_lst];
	end
    
    H = H(mask);
    Y = Y(mask);
    