%--------------------------------------------------------------------------
% varianload  - loads data files from our Varian spectrometer
%               amplitude corrected (sensitivity & nScans)
%
% Author: Anton Potocnik, F5, IJS
% Date:   25.01.2009
% Arguments:
%       [H, Y, exp] = varianload(filename)
%--------------------------------------------------------------------------

function [H, Y, exp] = varianload(filename) 

fid = fopen(filename);
    desc = textscan(fid, '%s %s',13,'headerLines', 1);  % read header, leave first row - troubles 
    fseek(fid, 0, 'bof');                               % goto the beginning of file
    C = textscan(fid, '%s %s','headerLines',14);        % read data
fclose(fid);

% Check for correct file format
if strcmp(desc{1,1},'MW_Frequency') == 0
    error(['Wrong file format: ' filename]);
end

% replace "," with "." and convert to double matrix
A = cell2mat(cellfun(@(x) str2double(strrep(x,',','.')),C,'UniformOutput',false));
H = A(:,1);
Y = A(:,2);

% replace "," with "." and convert to double for exp structure
freq = str2double(strrep(desc{1,2}{1},',','.'));  % MW_frequency
powr = str2double(strrep(desc{1,2}{2},',','.'));  % MW_power
cent = str2double(strrep(desc{1,2}{3},',','.'));  % Center_Field
swep = str2double(strrep(desc{1,2}{4},',','.'));  % Sweep
modb = str2double(strrep(desc{1,2}{5},',','.'));  % Modulation_Field
modf = str2double(strrep(desc{1,2}{6},',','.'));  % Modulation_Frequenfy
temp = str2double(strrep(desc{1,2}{7},',','.'));  % Temperature
ctim = str2double(strrep(desc{1,2}{8},',','.'));  % Conversation_Time
stim = str2double(strrep(desc{1,2}{9},',','.'));  % Sweep_Time
npnt = str2double(strrep(desc{1,2}{10},',','.')); % Number_Points
rcco = strrep(desc{1,2}{11},',','.');             % RC_Constant
sens = str2double(strrep(desc{1,2}{12},',','.')); % Sensitivity
nscn = str2double(strrep(desc{1,2}{13},',','.')); % Number_Scans

% exp.fname = filename;
exp.Temperature = temp;
exp.mwFreq = freq;
exp.Power = powr;
exp.CenterSweep = [cent swep]/10;   % easyspin units: [mT]
exp.Sensitivity = sens;
exp.nPoints = npnt;
exp.nScans = nscn;
exp.modFreq = modf;
exp.modB = modb;
exp.RC_Constant = rcco;
exp.cTime = ctim;

% Correct amplitude
Y = Y*sens/nscn;
    
    
    