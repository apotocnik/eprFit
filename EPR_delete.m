%--------------------------------------------------------------------------
% EPR_delete  - delete specters in epr structure
%
% Version: 1.0
% Author: Anton Potocnik, F5, IJS
% Date:   13.03.2009 - 15.03.2009
%       
% Arguments epr = EPR_delete(epr,v)
% Input:    
%       epr         
%       v           vector of indeces to delete
%--------------------------------------------------------------------------

function epr = EPR_delete(epr,v)

epr.data(v) = [];
epr.temp(v) = [];
epr.freq(v) = [];
if isfield(epr,'dates')
    epr.dates(v) = [];
end

if isfield(epr,'sim')
    if size(epr.sim,1) >= v
        epr.sim(v)=[];
    end
end
if isfield(epr.fit,'fits')
    if size(epr.fit.fits,1) >= v
        epr.fit.fits(v)=[];
    end
end
if isfield(epr.fit,'results')
    if size(epr.fit.results,1) >= v
        epr.fit.results(v,:)=[];
        epr.fit.results_g(v,:)=[];
    end
end
if isfield(epr.nra,'results')
    if size(epr.nra.results,1) >= v
    	epr.nra.results(v,:)=[];
        epr.nra.results_g(v,:)=[];
    end
end

epr.N = numel(epr.data);