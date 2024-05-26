%--------------------------------------------------------------------------
% EPR SIGNAL ANALYSIS Fit data using analitic functions
%
% Version: 2.0
% Author: Anton Potocnik @ IJS-F5
% Date:   27.10.2008 - 29.01.2009
% Input:  epr.fit.fitfun     ... string fit function
%         epr.fit.coef       ... structure with start values
%         epr.fit.options    ... parameters for fitting; start values will be
%                                updated
%         epr.fit.chosen     ... vector of specter numbers to fit
%         epr.fit.plot       ... plot results? (1/0)
%         epr.fit.range      ... specify range where to fit (extrange function)
% Output: epr.fit.fits{}.f           ... fit function
%                       .opts        ... fitting options
%                       .model       ... fitting model
%                       .results     ... 1.col coeff values, 2.col errors
%                       .gof         ... goodness of fit
%         epr.fit.results       ... all parameters with xc
%         epr.fit.results_g     ... all parameters with g (insted of xc)
%         epr.results_g#        ... parameters with g for a particular
%                                   funciton
%--------------------------------------------------------------------------

disp(' ');
disp('##########################################');
disp(['Fitting specters']);

clear Mov
%% Copy previous results
if isfield(epr.fit,'results')
    results = epr.fit.results;
end
if isfield(epr.fit,'results_g')
    results_g = epr.fit.results_g;
end

%% Get coefficient and problems names and values
coef_values = struct2cell(epr.fit.coef);
coef_names  = fieldnames(epr.fit.coef)';
param_names = coef_names;
coef_values = cell2mat(coef_values);

% Find vary and not vary indeces
vary_idx = find(coef_values(:,4)==1); 
not_vary_idx = find(coef_values(:,4)==0);  


% Extract not vary
prob_values = num2cell(coef_values(not_vary_idx,1))';
prob_names = coef_names(not_vary_idx);

% Extract vary
coef_values(not_vary_idx,:) = [];
coef_names(not_vary_idx) = [];

% Set fit model
epr.fit.model = fittype(epr.fit.fitfun,'problem',prob_names,'coefficients',coef_names);

epr.fit.options.StartPoint = coef_values(:,1)';        % Update start value
if strcmp(epr.fit.options.Algorithm,'Trust-Region')==1
    epr.fit.options.Upper = coef_values(:,3)';         % Update upper limits
    epr.fit.options.Lower = coef_values(:,2)';         % Update lower limits
else
    epr.fit.options.Upper = [];         % Delete upper limits
    epr.fit.options.Lower = [];         % Delete lower limits
end

%% FIT
model = epr.fit.model;
opts = epr.fit.options;

for i = epr.fit.chosen

    H = epr.data{i}.H;
    Y = epr.data{i}.Y;
    if isfield(epr.fit,'range')
        [H Y] = extrange(H,Y,epr.fit.range);
    end

    H = reshape(H,[],1);
    Y = reshape(Y,[],1);
    
    [f1 gof] = fit(H, Y, model, opts, 'problem', prob_values);
    coeff = coeffvalues(f1);        % Get fit coefficients
    
    %% Analyse and Save values
    epr.fit.fits{i}.f = f1;
    opts.StartPoint = coeff;
    epr.fit.fits{i}.opts = opts;
    epr.fit.fits{i}.model = model;
    epr.fit.fits{i}.gof = gof;
    
    % Extract errors
    ci = confint(f1,0.95);      % Boundaries within 95%
    err = (ci(2,:)-ci(1,:))/2;  % Absolut error
    
    % First column - coefficients
    res(vary_idx,1) = coeff;
    res(not_vary_idx,1) = cell2mat(prob_values);
    
    % Second column - errors
    res(vary_idx,2) = err;
    res(not_vary_idx,2) = zeros(1,numel(prob_values));

    epr.fit.fits{i}.results = res;
    epr.fit.fits{i}.coef = epr.fit.coef;

    
    %% Plot fit & data
    if epr.fit.plot==1
        figure(1);
        plot(H,Y,H,f1(H),'r');
        %axis([min(H) max(H) -0.01 0.01])
        legend(num2str(epr.temp(i)));
        Mov(i) = getframe;
    end
    
    tmp = sprintf('T=%3.0fK\tSSE=%3.6f\t', [epr.temp(i) gof.sse]);
%     disp(sprintf('%s%s', [tmp num2str(res(:,1))]));
    disp(sprintf('%s%s', [tmp num2str(res(:,1)')]));
    
end

%% Analyse Fits

for i = epr.fit.chosen
    vrstica = [epr.temp(i) reshape(epr.fit.fits{i}.results',1,[])];
    if exist('results','var')
        if size(results,2) ~= size(vrstica,2)   % different in size
            anse = questdlg('Different results already exist. Overwrite all previous results?','Overwrite');
            if strcmp(anse,'Yes')
                results = [];
                results_g = [];
            else
                return;
            end
        end 
    end
        
    results(i,:) = vrstica;
    results_g(i,:) = results(i,:);
    g_idx = 2*find(strncmpi('xc',param_names,2));  % 2* because of the errors
    for j = g_idx
        xc = results(i,j);
        dxc = results(i,j+1);
        g = xc2g(xc,epr.freq(i));
        results_g(i,j) = g;
        results_g(i,j+1) = g*dxc/xc;
    end
end

epr.fit.results = results;
epr.fit.results_g = results_g;
% epr.results_g = results_g;

%% Chop results to functions and plot individual function results
clear res

for i=1:50   % Search up to 5 functions
    r = strfind(param_names,num2str(i));
    r = cellfun('isempty',r);
    idx = find(r == 0);
    if isempty(idx); continue; end
    res = results_g(:,sort([1 2*idx 2*idx+1]));
%     epr.fit.(['results_g' num2str(i)]) = res;
    epr.(['results_g' num2str(i)]) = res;
    %% Plot
%     if epr.fit.plot==1
%         title_str = [epr.material '  ' epr.date '  ' num2str(i)];
%         plot_results(epr.fit.(['results_g' num2str(i)]),title_str,1);
%     end
end



clear i j H Y model opts f1 gof ind tmp coef_names coef_values prob_names prob_values ci err results coeff
clear vary_idx not_vary_idx  res param_names title_str xc dxc g g_idx results_g results r idx