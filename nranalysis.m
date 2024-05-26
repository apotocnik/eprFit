%--------------------------------------------------------------------------
% nranalysis  -  Numerical analysis of EPR spectra
%
% Author: Anton Potocnik, F5, IJS
% Date:   26.01.2009
% Arguments:
%       [Z, [A], [w], [xc]] = nranalysis(H, Y, [bline_corr], [int_cutoff],
%       [w_method], [xc_method])
% Input:
%        H, Y        - specter
%        bline_corr  - base line linear correction 
%                      [field_span_left field_span_right];
%                      the same units as in H
%        int_cutoff  - integrate only above int_cutoff (for calculating A)
%        w_method    - dH analysing method: 'p2p', 'fwhm', 'moment'
%        xc_method   - xc analysing method: 'p2p', 'fwhm', 'moment', 'max'
% Output:
%        Z  - specter integral
%        A  - Intensity (area under specter integral)
%        w  - dH
%        xc - xc
%--------------------------------------------------------------------------

%function [Z, A, w, xc] = nranalysis(H, Y, bline_corr, int_cutoff, w_method, xc_method)
function [Z, varargout] = nranalysis(H, Y, varargin)
bline_corr = 0;
int_cutoff = 0; 
w_method = 'p2p';
xc_method = 'p2p';
n = numel(varargin);
if n>0
    bline_corr = varargin{1};
end
if n>1
    int_cutoff = varargin{2};
end
if n>2
    w_method = varargin{3};
end
if n>3
    xc_method  = varargin{4};
end


%% First integral
H = reshape(H,[],1);
Z = cumtrapz(H,Y);

%% Integral Base Line Linear Correction
% Convert field values to index
if numel(bline_corr) == 2
   [C left_idx]  = min(abs(H-(H(1)+bline_corr(1))));
   [C right_idx] = min(abs(H-(H(end)-bline_corr(2))));
else
   [C left_idx]  = min(abs(H-(H(1)+bline_corr)));
   [C right_idx] = min(abs(H-(H(end)-bline_corr)));
end

if left_idx > 1 || right_idx < numel(H)     % Do nothing if bline_corr==0
    % Exctract data for linear fit
    x = H([1:left_idx right_idx:end]);
    z = Z([1:left_idx right_idx:end]);
    % Linear fit
    p = polyfit(x,z,1);

    % Correct base line
    Z = Z - polyval(p,H);
end

%% Integrate with cutoff
% if cutoff is empty dont cutoff
if isempty(int_cutoff)
    int_cutoff = -Inf;
end
A = trapz(H(Z>int_cutoff),Z(Z>int_cutoff));


%% Preanalysis
[ymin,iymin] = min(Y);
[ymax,iymax] = max(Y);
[zmax,izmax] = max(Z);


%% dH analysis
switch w_method
    case 'p2p'
        w = abs(H(iymin)-H(iymax));
    case 'fwhm'
        a = zmax/2;
        [C ileft] = min(abs(Z(1:izmax)-a));
        [C irigth] = min(abs(Z(izmax:end)-a));
        irigth = irigth + izmax;
        w = abs(H(irigth)-H(ileft));
    case 'moment'
        h = H(Z>int_cutoff);
        z = Z(Z>int_cutoff);
        xc = 1/A*trapz(h,z.*h);
        w = sqrt(1/A*trapz(h,z.*(h-xc).^2));
    otherwise
        error('Wrong analysis method! Try: p2p, fwhm, moment');
end


%% xc analysis
switch xc_method
    case 'p2p'
        xc = (H(iymin)+H(iymax))/2;
    case 'fwhm'
        if ~exist('ileft','var')
            a = zmax/2;
            [C ileft] = min(abs(Z(1:izmax)-a));
            [C irigth] = min(abs(Z(izmax:end)-a));
            irigth = irigth + izmax;
        end
        xc = (H(irigth)+H(ileft))/2;
    case 'moment'
        if ~exist('xc','var')      % If w_method != moment
            h = H(Z>int_cutoff);
            z = Z(Z>int_cutoff);
            xc = 1/A*trapz(h,z.*h);
        end
    case 'max'
        xc = H(izmax);
    otherwise
        error('Wrong analysis method! Try: p2p, fwhm, moment, max');
end

%% Output
varargout{1} = A;
varargout{2} = w;
varargout{3} = xc;

    
    