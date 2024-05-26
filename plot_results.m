%--------------------------------------------------------------------------
% plot_results  - plots multiple results: A(T), w(T), g(T)
%
% Author: Anton Potocnik, F5, IJS
% Date:   27.01.2009 - 15.03.2009
% Arguments:
%       h = plot_results(results_g,[title],[inv_A],[alpha])
% Input:    {results_g} ... [T, A, dA, w, dw, g, dg]   cell array of results
%           title       ... title string or epr structure
%           {alpha}     ... plot alpha ? (vector_alpha or 0) cell array
%--------------------------------------------------------------------------

function plot_results(results_g, title_str, alpha)
if nargin<2
    title_str = 'Title  Date';
end
if nargin<3
    alpha = 0;
end

if isstruct(title_str)
    if isfield(title_str,'material') && isfield(title_str,'date')
        title_str = [title_str.material '  ' title_str.date '  '];
    end
end

if ~iscell(results_g)
    results_g = {results_g};
end

if ~iscell(alpha)
    alpha = {alpha};
end

%--------------------------------------------------------------------------
font_size_title = 18;
font_size_labels = 14;
font_size_numbers = 14;
small_font_size_labels = 10;
small_font_size_numbers = 10;
%--------------------------------------------------------------------------

% Create figure
figure1 = figure('Position',[560,50,560,670]);

%% Create axes1 for A 
axes1 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'Position',[0.16 0.677 0.77 0.24],...
    'LineWidth',1.5,...
    'FontSize',font_size_numbers,...
    'FontName','Arial');
box('on');
hold('all');


% Create plot
for i=1:numel(results_g)
    % Get variables
    T = results_g{i}(:,1);
    A = results_g{i}(:,2);
    dA = results_g{i}(:,3);
    
    errorbar(T,A,dA,'Parent',axes1,'MarkerFaceColor',colors(i),...
    'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
end
grid on

% Create ylabel
ylabel('\chi_{EPR} (arb. units)','FontSize',font_size_labels,'FontName','Arial');

% Create title
title(title_str,'FontSize',font_size_title);

% % Create xlabel
% xlabel('\itT \rm(K)','FontSize',font_size_labels,'FontName','Arial');


%% Create axes2 for w
axes2 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'Position',[0.16 0.377 0.77 0.24],...
    'LineWidth',1.5,...
    'FontSize',font_size_numbers,...
    'FontName','Arial');
box('on');
hold('all');

% Create plot
for i=1:numel(results_g)
    % Get variables
    T = results_g{i}(:,1);
    w = results_g{i}(:,4);
    dw = results_g{i}(:,5);
    
    errorbar(T,w,dw,'Parent',axes2,'MarkerFaceColor',colors(i),...
    'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
end
grid on

% Create ylabel
ylabel('\Delta\itH\rm_{ } (mT)','FontSize',font_size_labels,'FontName','Arial');

% % Create xlabel
% xlabel('\itT \rm(K)','FontSize',font_size_labels,'FontName','Arial');


%% Create axes3 for g
axes3 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'Position',[0.16 0.077 0.77 0.24],...
    'LineWidth',1.5,...
    'FontSize',font_size_numbers,...
    'FontName','Arial');
box('on');
hold('all');

% Create plot
for i=1:numel(results_g)
    % Get variables
    T = results_g{i}(:,1);
    g = results_g{i}(:,6);
    dg = results_g{i}(:,7);
    
    errorbar(T,g,dg,'Parent',axes3,'MarkerFaceColor',colors(i),...
    'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
end
grid on

% Create ylabel
ylabel('\itg\rm-factor','FontSize',font_size_labels,'FontName','Arial');

% Create xlabel
xlabel('\itT \rm(K)','FontSize',font_size_labels,'FontName','Arial');


%% Create axes4 for 1/A
if 0
    axes4 = axes('Parent',figure1,'YMinorTick','off','XMinorTick','off',...
        'Position',[0.61 0.77 0.25 0.12],...
        'LineWidth',1.0,...
        'FontSize',small_font_size_numbers,...
        'FontName','Arial');
    box('on');
    hold('all');

    % Create plot
    for i=1:numel(results_g)
        % Get variables
        T = results_g{i}(:,1);
        A = results_g{i}(:,2);

        plot(T,1./A,'Parent',axes4,'MarkerFaceColor',colors(i),...
        'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
    end
    grid on

    % Create ylabel
    ylabel('1/\chi_{EPR} (a.u.)','FontSize',small_font_size_labels,'FontName','Arial');
    
    % Create xlabel
    xlabel('\itT \rm(K)','FontSize',small_font_size_labels,'FontName','Arial');
end


%% Create axes5 for alpha
if numel(alpha{1}) > 1
    axes5 = axes('Parent',figure1,'YMinorTick','off','XMinorTick','off',...
        'Position',[0.61 0.47 0.25 0.12],...
        'LineWidth',1.0,...
        'FontSize',small_font_size_numbers,...
        'FontName','Arial');
    box('on');
    hold('all');

    % Create plot
    for i=1:numel(alpha)
        T = results_g{i}(:,1);
        plot(T,alpha{i},'Parent',axes5,'MarkerFaceColor',colors(i),...
        'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
    end
    grid on

    % Create ylabel
    ylabel('alpha','FontSize',small_font_size_labels,'FontName','Arial');

    % Create xlabel
    xlabel('\itT \rm(K)','FontSize',small_font_size_labels,'FontName','Arial');
end
