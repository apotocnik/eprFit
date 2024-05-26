%--------------------------------------------------------------------------
% plot_fit_results  - plots multiple results from epr structre NEW!!!
%
% Author: Anton Potocnik, F5, IJS
% Date:   27.01.2009 - 13.10.2011
% Arguments:
%       h = plot_results(results_g,[title],[inv_A],[alpha])
%       h = plot_fit_results(epr)
%--------------------------------------------------------------------------

function plot_fit_results(epr)
% if nargin<2
%     title_str = 'Title  Date';
% end
% if nargin<3
%     alpha = 0;
% end

if ~isfield(epr,'fit')
    error('There is no epr.fit structure!')
end

if isfield(epr,'material') && isfield(epr,'date')
    title_str = [epr.material '  ' epr.date '  '];
end



% find parameters ---------------------------------------------------------
names = fieldnames(epr.fit.coef);

a={};
w={};
xc={};
alpha = {};
off=0;
for i=1:9
   txt = ['a' num2str(i)];
   in = find(strcmp(names,txt)==1);
   if in > 0 
       if i > 1
           off = in-1;
           in = in-off;
       end
       a{i} = 2*in;
   else
       break
   end;
   
   txt = ['w' num2str(i)];
   in = find(strcmp(names,txt)==1);
   if in > 0, w{i} = 2*(in-off); else w{i} = []; end;
   
   txt = ['xc' num2str(i)];
   in = find(strcmp(names,txt)==1);
   if in > 0, xc{i} = 2*(in-off); else xc{i} = []; end;
   
   txt = ['alpha' num2str(i)];
   in = find(strcmp(names,txt)==1);
   if in > 0, alpha{i} = 2*(in-off); else alpha{i} = []; end;
end

%--------------------------------------------------------------------------
% ask for additional plots

show_1chi = false;
show_alpha = false;

prompt = {'Show 1/chi?','Show alpha?','Show components:'};
dlg_title = 'Show additional plots';
num_lines = 1;
def = {'n','y',['1:' num2str(numel(a))]};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if strfind(answer{1},'y')
    show_1chi = true;
end
if strfind(answer{2},'y')
    show_alpha = true;
end
filter = str2num(answer{3});
rem = 1:numel(a);
rem(filter)=[];
if ~isempty(rem)
    a{rem} = [];
    w{rem} = [];
    xc{rem} = [];
    alpha{rem} = [];
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
for i=1:numel(a)
    % Get variables
    if isempty(a{i}), continue; end;
    T = epr.(['results_g' num2str(i)])(:,1);
    A = epr.(['results_g' num2str(i)])(:,a{i});
    dA = epr.(['results_g' num2str(i)])(:,a{i}+1);
    
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
for i=1:numel(w)
    % Get variables
    if isempty(w{i}), continue; end;
    T = epr.(['results_g' num2str(i)])(:,1);
    W = epr.(['results_g' num2str(i)])(:,w{i});
    dW = epr.(['results_g' num2str(i)])(:,w{i}+1);
    
    errorbar(T,W,dW,'Parent',axes2,'MarkerFaceColor',colors(i),...
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
for i=1:numel(xc)
    % Get variables
    if isempty(xc{i}), continue; end;
    T = epr.(['results_g' num2str(i)])(:,1);
    g = epr.(['results_g' num2str(i)])(:,xc{i});
    dg = epr.(['results_g' num2str(i)])(:,xc{i}+1);
    
    errorbar(T,g,dg,'Parent',axes3,'MarkerFaceColor',colors(i),...
    'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
end
grid on

% Create ylabel
ylabel('\itg\rm-factor','FontSize',font_size_labels,'FontName','Arial');

% Create xlabel
xlabel('\itT \rm(K)','FontSize',font_size_labels,'FontName','Arial');


%% Create axes4 for 1/A
if show_1chi
    axes4 = axes('Parent',figure1,'YMinorTick','off','XMinorTick','off',...
        'Position',[0.61 0.77 0.25 0.12],...
        'LineWidth',1.0,...
        'FontSize',small_font_size_numbers,...
        'FontName','Arial');
    box('on');
    hold('all');

    % Create plot
    for i=1:numel(a)
        % Get variables
        if isempty(a{i}), continue; end;
        T = epr.(['results_g' num2str(i)])(:,1);
        A = epr.(['results_g' num2str(i)])(:,a{i});


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
if numel(cell2mat(alpha)) > 0 && show_alpha
    axes5 = axes('Parent',figure1,'YMinorTick','off','XMinorTick','off',...
        'Position',[0.61 0.47 0.25 0.12],...
        'LineWidth',1.0,...
        'FontSize',small_font_size_numbers,...
        'FontName','Arial');
    box('on');
    hold('all');

    % Create plot
    for i=1:numel(alpha)
        if isempty(alpha{i}), continue; end;
        T = epr.(['results_g' num2str(i)])(:,1);
        alph = epr.(['results_g' num2str(i)])(:,alpha{i});
        plot(T,alph,'Parent',axes5,'MarkerFaceColor',colors(i),...
        'MarkerEdgeColor',colors(i),'Marker',markers(i),'LineStyle','none');
    end
    grid on

    % Create ylabel
    ylabel('alpha','FontSize',small_font_size_labels,'FontName','Arial');

    % Create xlabel
    xlabel('\itT \rm(K)','FontSize',small_font_size_labels,'FontName','Arial');
end


