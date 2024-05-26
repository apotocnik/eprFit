%--------------------------------------------------------------------------
% HansFormatload  - loads data files from HighField laboratory in Tallahassee
%
% Author: Anton Potocnik, F5, IJS
% Date:   25.01.2009-16.03.2010
% Arguments:
%       [H, X, Y, exp] = HansFormatload(filename)
%--------------------------------------------------------------------------

function [H, X, Y, exp] = HansFormatload(filename) 


% HEADER EXAMPLE:
% DataFile HansFormat 
% Fit of C:\Documents and Settings\Andrej\My Documents\LPS\Measurements\ESR\HFESR\Zn1Kagome\Zn1kagome_98pt4GHz_4Ku.dat column 1
% Fit to Lorentzians: 
% Offset : 0.00983316
% Linear term : -0.0114825
% Line          Amplitude    Center     Width    PhaseAngle(rad)  Splitting
% Line 1     0.0317099    3.21136    0.781157    -0.0188021    0.01
% Line 2     5.97433e-05    3.52008    0.0328154    -0.018    0.01
% a -0.00509997195427527    b 1796.57565175839
% >BEGIN

% Check for correct file format
desc = textread(filename,'%s',1);
if strcmp(desc{1,1},'DataFile') == 0
    error(['Wrong file format: ' filename]);
end

fid = fopen(filename);
    date = '';
    lowf = 0;
    higf = 0;
    sped = 0;
    temp = 0;
    sens = 1;
    modf = 0;
    modb = 0;
    rcco = 0;
    freq = 0;
    modb = 15; %Gauss : Typical value
    
    while ~strcmp(desc{1},'>BEGIN') % Skip the whole header - nothing useful
       desc = textscan(fid, '%s',1,'delimiter','\n');
    end

    C = textscan(fid, '%s %s %s %s %s');%,'headerLines',j);      % read data
fclose(fid);

% convert to double matrix
A = cell2mat(cellfun(@(x) str2double(x),C,'UniformOutput',false));
H = A(:,1);
X = A(:,2);
Y = A(:,4);
lowf = min(H);
higf = max(H);

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
    
    
    