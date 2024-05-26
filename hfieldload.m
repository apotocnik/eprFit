%--------------------------------------------------------------------------
% hFieldload  - loads data files from HighField laboratory in Tallahassee
%
% Author: Anton Potocnik, F5, IJS
% Date:   25.01.2009
% Arguments:
%       [H, X, Y, exp] = hfieldload(filename)
%--------------------------------------------------------------------------

function [H, X, Y, exp] = hfieldload(filename) 


% HEADER EXAMPLE:
% Measurement made Thu Jan 08 17:03:38 2009
% Main magnet sweep 11.80000 to 12.70000 Tesla at  2.001 mT/s 
% SweepCoil at   -0.000 mT 
% Main Magnet at 11.80000 T 
% Sweep Coil at -0.000000 T 
% ITC set at  300.000 K, T1: 300.000 K, T2: 509.100 K, T3: 1590.100 K 
% LockIn1: 500 uV, Mod Freq 50000.000000 Hz, Phase 54.700000 degr, Ampl. 5.000000 V
%   Tc 0.300000 s, Slope 18/oct, Dyn Res Normal .
% Source Freq 336 GHz, EIP lock 
% Vanadium Oxobenzoate
% 
% 300K
% 
% mod 16 Gauss
% >BEGIN

% Check for correct file format
desc = textread(filename,'%s',1);
if strcmp(desc{1,1},'Measurement') == 0
    error(['Wrong file format: ' filename]);
end

fid = fopen(filename);
    % read date
%     desc = textscan(fid, '%s','delimiter','\n');
    
    desc = textscan(fid, 'Measurement made %s %s %s %s %s\n');   
    date = '';
    date = [desc{1}{1} ' ' desc{2}{1} ' ' desc{3}{1} ' ' desc{4}{1} ' ' desc{5}{1}];
    
    % read scan range
    desc = textscan(fid, 'Main magnet sweep %s to %s Tesla at %s mT/s');
    
    lowf = str2double(desc{1});
    higf = str2double(desc{2});
    sped = str2double(desc{3});
    
    desc = textscan(fid, '%s',3,'delimiter','\n');      % Skip 3 lines
    
    desc = textscan(fid, 'ITC set at %f K, T1: %f K, T2: %f K, T3: %f K \n');
    temp = desc{2};
    
    desc = textscan(fid, 'LockIn1: %f %s Mod Freq %f Hz, Phase %f degr, Ampl. %f V\n');
    sens = desc{1}*units2value(desc{2}{1}(1));
            
    modf = desc{3}/1000;  % units [kHz]
    modb = desc{5};
    
    %desc = textscan(fid, '  Tc %f s, Slope %f/oct, Dyn Res Normal .\n'); %
    % can have some other text
    % skip 2x \n
    desc = textscan(fid, '%s',2,'delimiter','\n');
    rcco = 0; % sscanf ... desc{1};
    
    while ~strcmp(desc{1,1},'Source') % Skip all the words and come to Source Freq ...
       desc = textscan(fid, '%s',1);
    end
    
    desc = textscan(fid, 'Freq %f GHz, EIP lock \n');
    freq = desc{1};
    if isempty(freq)
        answer = inputdlg('Frequency ? (GHz)','Error reading freqeuncy',1,{'0'});
        freq=str2double(answer);
    end
    
    modb = 15; %Gauss : Typical value
    
    while ~strcmp(desc{1},'>BEGIN')
       desc = textscan(fid, '%s',1,'delimiter','\n');
    end

    C = textscan(fid, '%s %s %s');%,'headerLines',j);      % read data
fclose(fid);

% convert to double matrix
A = cell2mat(cellfun(@(x) str2double(x),C,'UniformOutput',false));
H = A(:,1);
X = A(:,2);
Y = A(:,3);

% exp.fname = filename;
exp.Temperature = temp;
exp.mwFreq = freq;
exp.CenterSweep = [(lowf+higf)/2 (higf-lowf)]*1000;   % easyspin units: [mT]
exp.Sensitivity = sens;
exp.nPoints = numel(H);
%exp.nScans = 1;
exp.modFreq = modf;     % units: [kHz]
exp.modB = modb;        % units: [G]
exp.RC_Constant = rcco;
exp.cTime = exp.CenterSweep(2)/sped;

% No need to correct amplitude; allways in mV !!!
% Y = Y*sens;
% X = X*sens;
    
    
    