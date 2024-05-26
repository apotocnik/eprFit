%--------------------------------------------------------------------------
% epr_sort  - sorts specters in epr structure
%
% Version: 1.1
% Author: Anton Potocnik, F5, IJS
% Date:   17.05.2009
%       
% Input:    
%       epr         
%       epr.sort    temp,freq,dates
%       epr.sortDir descend,ascend
%--------------------------------------------------------------------------


if ~isfield(epr,'sortDir')
    epr.sortDir = 'descend';  % 'ascend'
end

% SORT in ascend order
[v ix] = sort(epr.(epr.sort),epr.sortDir);


% if strcmp(epr.sortDir,'descend')    % problems with cell arrays
%     ix = ix(numel(ix):-1:1);
% end

epr.data = reshape(epr.data(ix),[],1);
epr.temp = reshape(epr.temp(ix),[],1);
epr.freq = reshape(epr.freq(ix),[],1);
epr.dates = reshape(epr.dates(ix),[],1);

epr.fit.fits = reshape(epr.fit.fits(ix),[],1);
epr.sim = reshape(epr.sim(ix),[],1);

if isfield(epr.fit,'results')
    s = size(epr.fit.results,1);
    if s==0, return; end
    epr.fit.results(s+1:epr.N,:) = zeros(epr.N-s,size(epr.fit.results,2));
    epr.fit.results = epr.fit.results(ix,:);
end

if isfield(epr.fit,'results_g')
    s = size(epr.fit.results_g,1);
    if s==0, return; end
    epr.fit.results_g(s+1:epr.N,:) = zeros(epr.N-s,size(epr.fit.results_g,2));
    epr.fit.results_g = epr.fit.results_g(ix,:);
end

if isfield(epr,'results_g1')
    s = size(epr.results_g1,1);
    if s==0, return; end
    epr.results_g1(s+1:epr.N,:) = zeros(epr.N-s,size(epr.results_g1,2));
    epr.results_g1 = epr.results_g1(ix,:);
end

if isfield(epr,'results_g2')
    s = size(epr.results_g2,1);
    if s==0, return; end
    epr.results_g2(s+1:epr.N,:) = zeros(epr.N-s,size(epr.results_g2,2));
    epr.results_g2 = epr.results_g2(ix,:);
end

if isfield(epr,'results_g3')
    s = size(epr.results_g3,1);
    if s==0, return; end
    epr.results_g3(s+1:epr.N,:) = zeros(epr.N-s,size(epr.results_g3,2));
    epr.results_g3 = epr.results_g3(ix,:);
end

if isfield(epr,'results_g4')
    s = size(epr.results_g4,1);
    if s==0, return; end
    epr.results_g4(s+1:epr.N,:) = zeros(epr.N-s,size(epr.results_g4,2));
    epr.results_g4 = epr.results_g4(ix,:);
end

if isfield(epr.nra,'results')
    s = size(epr.nra.results,1);
    if s==0
        return
    end
    if s < epr.N
       epr.nra.results(s:epr.N,:)=zeros(size(epr.nra.results,2),epr.N-s);
       epr.nra.results_g(s:epr.N,:)=zeros(size(epr.nra.results_g,2),epr.N-s);
    end
    epr.nra.results = epr.nra.results(ix,:);
    epr.nra.results_g = epr.nra.results_g(ix,:);
end

clear v ix
    