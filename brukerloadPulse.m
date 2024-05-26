%--------------------------------------------------------------------------
% brukerload  - loads data files from our Bruker spectrometer
%               amplitude corrected (sensitivity & nScans)
%
% Author: Anton Potocnik, F5, IJS
% Date:   12.03.2011
% Arguments:
%       [H, Y, exp] = brukerload(filename)
%--------------------------------------------------------------------------

function [H, Y, exp] = brukerloadPulse(filename) 

    [Path Name Ext] = fileparts(filename);
    Path = [Path '\'];
    
    % Read .DSC description file
    id = fopen([Path Name '.DSC']);
    if id == -1, error([Name ' cannot be found !!!']); end;
    str = fgetl(id);
    while ischar(str)
        if strfind(str,'DATE');
            DATE = str(9:end);
        end
        if strfind(str,'TIME');
            TIME = str(9:end);
        end
        if strfind(str,'MWFQ');             % MW_frequency
            freq = str2double(str(9:end));
        end
        if strfind(str,'MWPW');             % MW_power
            powr = str2double(str(9:end));
        end
        if strfind(str,'CenterField');      % Center Field
            cent = str2double(str(15:end-2));
        end
        if strfind(str,'SweepWidth');       % Sweep
            swep = str2double(str(14:end-2));
        end
        if strfind(str,'ModAmp');           % Modulation_Field
            modb = str2double(str(10:end-2));
        end
        if strfind(str,'ModFreq');           % Modulation_Field
            modf = str2double(str(11:end-4));
        end
        if strfind(str,'TimeConst');         % RC_Constant
            rcco = str2double(str(15:end-3));
        end
        if strfind(str,'SweepTime');         % Sweep_Time
            stim = str2double(str(11:end-2));
        end
        if strfind(str,'NbScansDone');       % Number of scans
            nscn = str2double(str(15:end));
        end
        if strfind(str,'Resolution ');        % Number of points
            npnt = str2double(str(14:end));
        end
        if strfind(str,'RCAG');              % Sensitivity
            RCAG = str2double(str(5:end));
        end
        if strfind(str,'ConvTime');              % Sensitivity
            ctim = 1e-3*str2double(str(14:end-3));
        end
         
       str = fgetl(id);
    end
    fclose(id);

    % Read .DSC description file
    id = fopen([Path Name '.DTA'],'r','ieee-be');
    if id == -1, error([name ' cannot be found !!!']); end;
    Y = fread(id,[npnt,1],'float64','ieee-be');
    fclose(id);
    
%     brukerread([Path Name '.DTA'])

    H = linspace((cent-swep/2),(cent+swep/2),npnt); % in Gauss

    % exp.fname = filename;
    % exp.Temperature = temp;
    exp.mwFreq = freq;
    exp.Power = powr;
    exp.CenterSweep = [cent swep]/10;   % units: [mT]
%     exp.Sensitivity = sens;
    exp.nPoints = npnt;
    exp.nScans = nscn;
    exp.modFreq = modf;
    exp.modB = modb;
    exp.RC_Constant = rcco;
    exp.cTime = ctim;
    exp.RCAG = RCAG;
    
    % Correct amplitude
    Y = Y/nscn;
    
    
    