%--------------------------------------------------------------------------
%% EPR SIGNAL ANALYSIS Load data
% Version 2.1 - Loads almost everything
%
% Author: Anton Potocnik, F5, IJS
% Date:   27.10.2008 - 22.03.2010
% Input:  epr.path 
%         [epr.sort]    % 'temp' 'date'
% Output  epr.data{}.exp .fname, .Temperature, .Frequency, .Power,
%                        .CenterSweep, .Sensitivity, .nPoints, .nScans, 
%                        .modFreq, .modB, .RC_Constant, .cTime
%                   .H   (magnetic field (mT)) !!!
%                   .Y   (dP/dH (arb. units))
%                   .X 
%            .temp()
%            .dates()
%            .freq()
%            .N         % numel(epr.data)
%--------------------------------------------------------------------------

%% Get file list and its information
files = dir(epr.path);
fpath = fileparts(epr.path);

if isempty(files)
    error('Files not found!')
end

%% Get previous data
ind = 1;   % data counter, different than i if data allready exist
if isfield(epr,'data')      % if data allready exists add new data
    ind = numel(epr.data)+1;
    data = epr.data;
    temp = epr.temp;
    date = epr.dates;
    freq = epr.freq;
end

%% Load files
disp(' ');
disp('##########################################');
disp('Loading data from:');
disp(['  ' fpath]);

for i=1:numel(files)
    if files(i).isdir == true  % Skip Folders
        continue;
    end
    fname = fullfile(fpath,files(i).name);
    
    % Recognize file format
    desc = textread(fname,'%s',1,'bufsize',16380); 
    
%     if ~isempty(str2num(desc{1})) % Plain column ascii numbers
%         Ma = dlmread(fname);
%         if size(Ma,2) == 2
%             H = Ma(:,1);
%             Y = Ma(:,2);
%             X = zeros(size(Y));
%         elseif size(Ma,2) == 3
%             H = Ma(:,1);
%             Y = Ma(:,2);
%             X = Ma(:,3);
%         end
%         clear Ma
%         exp.mwFreq = 0;
%         exp.Temperature = 295;
%         exp.nPoints = numel(H);
%         exp.CenterSweep = [mean(H) (H(end)-H(1))];
%     else
        switch desc{1}
            case 'Sample_Description:'              
                [H Y exp] = varianload(fname);      % Varian
                H = H/10;              % Units are now in mT
                if isfield(exp,'CenterField')
                    exp.CenterField = exp.CenterField/10; % mT
                end
            case 'Measurement'                      
                [H X Y exp] = hfieldload(fname);    % HighField Tallahassee
                H = H*1e3;             % Units now in mT
            case 'DataFile'
                [H X Y exp] = HansFormatload(fname);    % HighField Tallahassee
                H = H*1e3;             % Units now in mT
            case 'SIMPLE'
                [Y,H] = ReadFits(fname);
                exp.mwFreq = 0;
                exp.Temperature = 295;
            case '#DESC'
                [H Y expe] = brukerload(fname);         % EasySpin: Bruker, ...
                H = H/10;              % Units are now in mT
                expe.mwFreq = expe.mwFreq/1e9;          % Units: GHz
                if isfield(expe,'CenterField')
                    expe.CenterSweep = [str2double(strtok(expe.CenterField,' G')) str2double(strtok(expe.SweepWidth,' G'))];
                    expe.CenterField = expe.CenterField/10; % mT
                end
                expe.Temperature = 300;
                if isfield(expe,'RCAG')
                    expe.RCAG = 50; % Added bz AP to awoid confution. All this part should be removed!!!
                    Y = Y*exp(-expe.RCAG/20*log(10));  % function exp has the same name as variable exp!
                end
                exp=expe;
                
            otherwise
                [Ma, delimiter, nheaderlines] = importdata(fname);
                prompt = ['Which columns to use for H, Y, [X] ?\n\n'];
                for j = 1:nheaderlines
                    prompt = [prompt Ma.textdata{j} '\n'];
                end
                prompt = [prompt '\nH:'];
                answer = inputdlg({sprintf(prompt),'Y:','X:'},'Import data',1,{'1','2',''});
%                 Ma = dlmread(fname);

                if isempty(answer), return; end
                Ma = Ma.data;
                if isempty(answer{3})
                    H = Ma(:,str2double(answer{1}));
                    Y = Ma(:,str2double(answer{2}));
                    X = zeros(size(Y));
                else
                    H = Ma(:,str2double(answer{1}));
                    Y = Ma(:,str2double(answer{2}));
                    X = Ma(:,str2double(answer{3}));
                end
                clear Ma
                exp.mwFreq = 0;
                exp.Temperature = 295;
                exp.nPoints = numel(H);
                exp.CenterSweep = [mean(H) (H(end)-H(1))];
                
        end
%     end
    
    H = reshape(H,[],1);
    Y = reshape(Y,[],1);
    
%     data{ind}.fname = files(i).name;
    tmp = pwd;
    cd(fpath);
    data{ind}.fname = fullfile(pwd,files(i).name);
    cd(tmp);
    clear tmp
    data{ind}.date = files(i).date;
    
    if ~isreal(Y)
        X = imag(Y);
        Y = real(Y);
    end
    
    if ~exist('X','var')  % Create X if it does not exist!
        X = zeros(size(Y));
    end
    
    pos = find(isnan(H));
    if ~isempty(pos)
       H(pos) = [];
       Y(pos) = [];
       X(pos) = [];
    end
    
    % Sort field to allways increasing
    [v ix] = sort(H);   
    data{ind}.H = H(ix);
    data{ind}.Y = Y(ix);
    data{ind}.X = X(ix);
    if isfield(exp,'CenterSweep')
        exp.CenterSweep = abs(exp.CenterSweep);
    end
    
    data{ind}.exp = exp;
    
    if isfield(exp,'Temperature')
        temp(ind) = exp.Temperature;
    else
        temp(ind) = 295;     % Room temperature 22°C
    end
    
    if isempty(exp.mwFreq)
        freq(ind) = 0;
    else
        freq(ind) = exp.mwFreq;
    end
    
    if isfield(exp,'DATE')
        date{ind} = [exp.DATE ' ' exp.TIME];
    else
        date{ind} = files(i).date;
    end
    
    
    
    disp(sprintf('%d\t%s\tT=%3.2fK\t f=%3.2fGHz',ind,files(i).name,temp(ind),freq(ind)));
    ind = ind+1;
end

%% SORT and write to epr structure
% if isfield(epr,'sort')
%     switch epr.sort
%         case 'dates'
%             [v ix] = sort(dates);
%         case 'freq'
%             [v ix] = sort(freq);
%         otherwise
%             [v ix] = sort(temp,'ascend');
%     end
% else 
%     [v ix] = sort(temp,'descend');
% end
% epr.data = reshape(data(ix),[],1);
% epr.temp = reshape(temp(ix),[],1);
% epr.freq = reshape(freq(ix),[],1);
% epr.dates = reshape(date(ix),[],1);

epr.data = data;
epr.temp = temp;
epr.freq = freq;
epr.dates = date;

epr.N = numel(epr.data);

clear temp freq date X Y H i ind files fpath v ix data desc exp fname

