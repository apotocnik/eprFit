% Function for reading and plotting data from fits files generated by Andor
% Shamrock spectrometer and Andor Newton EMCCD camera.
%
%
% USAGE:
% [Data,Wavelength,FileInfo]=ReadFits(filename);
% [Data,Wavelength,FileInfo]=ReadFits(filename,'PropertyName',PropertyValue,...);
%
% INPUT arguments:
%  - filename: Filename or just a part of filname of fits file.
%  - optional property name/value pairs (see below)
% PROPERTIES:
%  - xAxis3D = AqNo | {Time} | CustTime | CustStep
%    In what units is the x axis in kinetic series.
%       - AqNo: Displays aquisition number.
%       - Time: Original time, if not avaiable uses AqNo.
%       - CustTime: Custom units calculated from time. Requires
%                   xAxisStart, xAxisEnd, elapsedTime and xAxisLabel.
%       - CustStep: Custom units calculated from aquisition number.
%                   Requires xAxisStart, xAxisStep and xAxisLabel.
%  - xAxisStart = double
%    Initial value of the changing parameter (e.g. temperature, voltage)
%  - xAxisEnd = double
%    Final value of the changing parameter (e.g. temperature, voltage)
%  - elapsedTime = double
%    Total elapsed time in kinetic series in seconds. Needed to calculate
%    custom x axis.
%  - xAxisStep = double
%    What is the step/change of the changing parameter in one aquisition.
%  - xAxisLabel = string
%    Label of the custom axis.
%  - ROI = vector of length=2, {[1:200]}
%    In the case of kinetic series of images or just one image the ROI 
%    defines the horizontal stripe where to integrate to get spectrum.
%  - SingleSpectrum = double
%    In the case of kinetic aquisition (either FVB or Image) lets you
%    display just one aquisiton by aquisition number.
%  - BaselineCorr = 0 | {1}
%    Whether to substract the background.
%    Default value for single spectrum is 0.
%  - NonLinInt = {linear}, sqrt, log, NumericValue
%    How is intensity displayed. In the case of NumericValue
%    Data=Data.^NonLinInt; For 3D the default is 0.8.
%  - CutLowInt = {0} | 1
%    Whether to cut off aquisitions with too low intensity.
%
% OUTPUT arguments:
%  - Data: Matrix containing modified measurment data.
%  - Wavelength: Vector containing wavelength data.
%  - FileInfo: Information about the file generated by fitsinfo.
%
% EXAMPLE:
% ReadFits('s05','xAxis3D','CustTime','xAxisStart',300,'xAxisEnd',340,'elapsedTime',40*60,'xAxisLabel','Temperature (\circC)','CutLowInt',1,'NonLinInt','log');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Matjaz Humar                                 %
%                       Jozef Stefan Institute                            %
%                         matjaz.humar@ijs.si                             %
%                            Version 3.1.3 - LAST VERSION                 %
%                                8.8.2010                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Data,wavelength,FileInfo]=ReadFits(filename,varargin)

%--------------------------------------------------------------------------
% READ PROPERTY NAME/VALUE PAIRS

% xAxis3D: AqNo, Time, CustTime, CustStep
% Default is Time or if no time info available then AqNo
xAxis3D=PropertyValue(varargin,'xAxis3D','Default');

% For custom x axis based on time we need 4 parameters
if strcmp(xAxis3D,'CustTime')
    
    xAxisStart=PropertyValue(varargin,'xAxisStart','PropertyRequired');
    xAxisEndRelative=PropertyValue(varargin,'xAxisEnd','PropertyRequired');
    elapsedTime=PropertyValue(varargin,'elapsedTime','PropertyRequired');
    xAxisLabel=PropertyValue(varargin,'xAxisLabel','PropertyRequired');
    
end


% For custom x axis based on aqusition number we need 3 parameters
if strcmp(xAxis3D,'CustStep')
    
    xAxisStart=PropertyValue(varargin,'xAxisStart','PropertyRequired');
    xAxisStep=PropertyValue(varargin,'xAxisStep','PropertyRequired');
    xAxisLabel=PropertyValue(varargin,'xAxisLabel','PropertyRequired');
    
end

% Used in image/kinetics to define region of interest
% ROI: [startRow,endRow]
% Default [1,200] - full verical binning (for image/kintetics)
% For single image it displays it as default as an image
ROI=PropertyValue(varargin,'ROI','Default');


% Used in kinetics to define just one aquisition
% SingleSpectrum: aqusition number
% Default 0 - whole kinetic series
SingleSpectrum=PropertyValue(varargin,'SingleSpectrum',0);


% BaselineCorr: 0 or 1
% Default 1
BaselineCorr=PropertyValue(varargin,'BaselineCorr','Default');


% NonLinInt: linear, sqrt, log
% Default: linear for single spectrum, 0.8 for 3D
NonLinInt=PropertyValue(varargin,'NonLinInt','Default');


% Cut low intensity aquisitions
% CutLowInt: 0 or 1
% Default 0
CutLowInt=PropertyValue(varargin,'CutLowInt',0);


% READ PROPERTY NAME/VALUE PAIRS - END
%--------------------------------------------------------------------------


% Contruct full filename from partial filename
%filenameFull=readFilename(filename);
filenameFull = filename;

% Read data and file info
Data = (fitsread(filenameFull));
FileInfo=fitsinfo(filenameFull);

% Data containing various information
Keywords=strvcat(FileInfo.PrimaryData.Keywords{:,1}); %#ok<VCAT>

% Find and read the calibration data
CalibPosition = strmatch('CALIB', Keywords);
Calibration = str2double(strread(FileInfo.PrimaryData.Keywords{CalibPosition,2}, '%s', 'delimiter', ','));


% Find and read the acquisition mode ('Single Scan' or 'Accumulate' or 'Kinetics')
AcquisitionModePosition = strmatch('ACQMODE', Keywords);
AcquisitionMode = strread(FileInfo.PrimaryData.Keywords{AcquisitionModePosition,2}, '%s', 'delimiter', ',');
% AcquisitionMode='Single Scan';

% Find and read the readout mode ('Full Vertical Binning' or 'Image')
ReadoutModePosition = strmatch('READMODE', Keywords);
ReadoutMode = strread(FileInfo.PrimaryData.Keywords{ReadoutModePosition,2}, '%s', 'delimiter', ',');

%--------------------------------------------------------------------------

% Interchange first two dimensions of the matrix
% Final result examples:
%    - FVB/single: 1600x1
%    - FVB/kinetic: 1600x10
%    - image/single: 1600x200
%    - image/kinetic: 1600x200x10
Data=permute(Data,[2 1 3]);


% In the case we want only one aquisition from kinetic series, this
% transforms the series to just one spectrum or image
if SingleSpectrum>0
    if strcmp(ReadoutMode,'Full Vertical Binning')
        Data=Data(:,SingleSpectrum);
    elseif strcmp(ReadoutMode,'Image')
        Data=Data(:,:,SingleSpectrum);
    end
    
    AcquisitionMode='Single Scan';
end


% In the case we want a single spectra from single image it uses ROI for vertical integration
if strcmp(ReadoutMode,'Image')
    if ~strcmp(ROI,'Default')
        Data=sum(Data(:,ROI(1):ROI(2)),2);
        ReadoutMode='Full Vertical Binning';
    end
end



% Find out the wavelength and kinetic axis lenght
DataSize=FileInfo.PrimaryData.Size;
xAxisLength=DataSize(1,1); % First dimension is wavelength length
YAxisLength=DataSize(1,ndims(Data)); % Last dimension is kinetic axis lenght



% Calculate wavelength from calibration
wavelength=...
    Calibration(1).*ones(xAxisLength,1)+...
    Calibration(2).*[1:xAxisLength]'   +...
    Calibration(3).*[1:xAxisLength]'.^2+...
    Calibration(4).*[1:xAxisLength]'.^3;    %#ok<NBRAK>

%--------------------------------------------------------------------------

% Sets the xAxis in kinetics acquisition mode
if strcmp(xAxis3D,'AqNo')
    xAxis=1:YAxisLength;
    xAxisLabel='Acquisition';
    % Find and read the time step
elseif strcmp(xAxis3D,'Time') || strcmp(xAxis3D,'CustTime') || strcmp(xAxis3D,'Default')
    TimeStepPosition = strmatch('KCT', Keywords);
    if ~isempty(TimeStepPosition)
        TimeStep = (FileInfo.PrimaryData.Keywords{TimeStepPosition,2});
        % Calculate time from time step
        if strcmp(xAxis3D,'Time') || strcmp(xAxis3D,'Default')
            xAxis=0:TimeStep:TimeStep*(YAxisLength-1);
            xAxisLabel='Time (s)';
        elseif strcmp(xAxis3D,'CustTime')
            xAxisStep=(xAxisEndRelative-xAxisStart)/elapsedTime*TimeStep;
            xAxisEnd=xAxisStart+xAxisStep*(YAxisLength-1);
            xAxis=xAxisStart:xAxisStep:xAxisEnd;
        end
    else
        if strcmp(xAxis3D,'Time') || strcmp(xAxis3D,'CustTime')
            disp('Warning: Time data was not found. Aquisition number will be used instead.')
        end
        xAxis=1:YAxisLength;
        xAxisLabel='Acquisition';
    end
elseif strcmp(xAxis3D,'CustStep')
    xAxis=xAxisStart:xAxisStep:YAxisLength*xAxisStep;
end

%--------------------------------------------------------------------------



if strcmp(AcquisitionMode,'Single Scan') || strcmp(AcquisitionMode,'Accumulate')
    
    if strcmp(ReadoutMode,'Full Vertical Binning')
        if strcmp(BaselineCorr,'Default')
            BaselineCorr=0;
        end
        Data=BaselineCorrFunc(Data,BaselineCorr);
        Data=NonLinIntFunc(Data,NonLinInt);
        plot(wavelength,Data)
        xlabel('Wavelength (nm)','FontSize',16)
        ylabel('Intensity (Arbitrary Units)','FontSize',16)
    elseif strcmp(ReadoutMode,'Image')
        Data=NonLinIntFunc(Data,NonLinInt);
        %scrsz = get(0,'ScreenSize');
        figure('Position',[50 500 1150 300])
        surf(wavelength,1:YAxisLength,Data')
        view([0 90]);
        xlim([min(wavelength),max(wavelength)])
        xlabel('Wavelength (nm)','FontSize',16)
        shading flat
    else
        error(['Unrecognized readout mode: ' ReadoutMode])
    end
    
elseif strcmp(AcquisitionMode,'Kinetics')
    
    if strcmp(ReadoutMode,'Full Vertical Binning')
        if strcmp(BaselineCorr,'Default')
            BaselineCorr=1;
        end
        Spectra3D(wavelength,xAxis,xAxisLabel,Data,BaselineCorr,NonLinInt,CutLowInt)
    elseif strcmp(ReadoutMode,'Image')
        % Sum in vertical direction
        if strcmp(ROI,'Default')
            ROI=[1,200];
        end
        Data2=sum(Data(:,ROI(1):ROI(2),:),2);
        Data=[];
        Data(:,:)=Data2(:,:,:);
        if strcmp(BaselineCorr,'Default')
            BaselineCorr=1;
        end
        Spectra3D(wavelength,xAxis,xAxisLabel,Data,BaselineCorr,NonLinInt,CutLowInt)

    else
        error(['Unrecognized readout mode: ' ReadoutMode])
    end
    
else
    
    error(['Unrecognized acquisition mode: ' AcquisitionMode])
    
end

% Set larger font
set(gca,'FontSize',14)


end

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


function Value=PropertyValue(Arguments,PropertyName,DefaultValue)

[PropertyFound,PropertyPos]=max(strcmp(Arguments,PropertyName));
if PropertyFound==1
    Value=Arguments{PropertyPos+1};
else
    % If the property is missing either put the default value or
    % if the property is required return error.
    if strcmp(DefaultValue,'PropertyRequired')
        error(['Parameter ' PropertyName ' missing.'])
    else
        Value=DefaultValue;
    end
end

end


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


% Recognizes filename with or without extension. Filename can also include
% just couple of unique characters.
function filenameFull=readFilename(filename)

filenameFull=[];
% All the .fits files
AllFiles=dir('*.fits');

for k=1:size(AllFiles,1)
    
    CurrentFilename=AllFiles(k).name;
    if strfind(CurrentFilename,filename)
        % If a filename was already found, then more files with the same string exist.
        if isempty(filenameFull)
            filenameFull=CurrentFilename;
        else
            error('More files with the same string exist. Type in full filename.')
        end
    end
end
end


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


function Spectra3D(wavelength,xAxis,xAxisLabel,Data,BaselineCorr,NonLinInt,CutLowInt)

% Function for drawing kinetic series of spectra
% USAGE: Spectra3D(wavelength,xAxis,xAxisLabel,Data,baseline)
% INPUT arguments:
%  - wavelength: Vector containing wavelength data.
%  - xAxis: Vector cointaining x axsis data.
%  - xAxisLabel: x axis label.
%  - Data: Martrix containing raw data from fits file
%  - baseline: Weather the background should be substracted.
% OUTPUT arguments: NONE
% EXAMPLE: Spectra3D(wavelength,xAxis,xAxisLabel,Data,1)


% Cut off the part with zero intensity (if the kinetic aquisition was
% prematurely stopped) or where the intensity is lower than certain
% treshold
for i=1:size(Data,2)
    if sum(Data(:,i))<10000 && CutLowInt
        disp(['Warning: The kinetic series was cut off at ' num2str(i)...
            ' of the total ' num2str(size(Data,2)) ' acquisitions, due to too low intensity.'])
        Data(:,i:end)=[];
        xAxis(i:end)=[];
        break
    end
end


% Cut off part of the kinetic series because of memory limitations.
maxKinteicSize=3000;
if size(Data,2)>maxKinteicSize
    disp(['Warning: The kinetic series was cut off at ' num2str(maxKinteicSize)...
        ' of the total ' num2str(size(Data,2))...
        ' acquisitions. Max allowed kinetic series size is ' num2str(maxKinteicSize) '.'])
    Data(:,maxKinteicSize+1:end)=[];
    xAxis(maxKinteicSize+1:end)=[];
end

Data=BaselineCorrFunc(Data,BaselineCorr);

% For 3D the default nonlinarity is 0.8
if strcmp(NonLinInt,'Default')
    NonLinInt=0.8;
end
Data=NonLinIntFunc(Data,NonLinInt);

% Additionally smooth the data
%Data=smooth2a(Data,1,1);

fullscreen = get(0,'ScreenSize');
figure('Position',[5 35 fullscreen(3)-35 fullscreen(4)-110])
pcolor(xAxis,wavelength,Data);
xlabel(xAxisLabel,'FontSize',16)
ylabel('Wavelength (nm)','FontSize',16)
shading flat
%colormap hot

end

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

function Data=BaselineCorrFunc(Data,BaselineCorr)

% Substract the background
if BaselineCorr
    % In case of single spectrum use 1D smoothing
    if size(Data,2)==1
        background=smooth(Data,1000,'rloess');
        Data=Data-background;
        Data=Data-min(Data);
    % In case of kinetic aquisition use 2D smoothing
    elseif size(Data,2)>1
        background=smooth2a(Data,300,300);
        Data=Data-background;
        Data=Data-min(min(Data));
    else
        disp('Something wrong with the size of the data.')
    end
end

end

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

function Data=NonLinIntFunc(Data,NonLinInt)
% Nonlinear intensity scale
if isa(NonLinInt,'numeric')
    Data=Data.^NonLinInt;
elseif strcmp(NonLinInt,'linear') || strcmp(NonLinInt,'Default')
    % Do not do anything.
elseif strcmp(NonLinInt,'sqrt')
    Data=sqrt(Data);
elseif strcmp(NonLinInt,'log')
    Data=log(Data);
end

end



%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


function matrixOut = smooth2a(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc)
%
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element "i" by the mean of the rectange centered on "i".  Any NaN
% elements are ignored in the averaging.  If element "i" is a NaN, then it
% will be preserved as NaN in the output.  At the edges of the matrix,
% where you cannot build a full rectangle, as much of the rectangle that
% fits on your matrix is used (similar to the default on Matlab's builtin
% function "smooth").
%
% "matrixIn": original matrix
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
%
% "matrixOut": smoothed version of original matrix
%
%
% 	Written by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
%
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
%
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se

%
% Initial error statements and definitions
%
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr;
if nargin < 3, N(2) = N(1); else N(2) = Nc; end

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the
% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by-
% (2*Nc+1) rectangle centered on element "i".
%
[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

%
% Setting all "NaN" elements of "matrixIn" to zero so that these will not
% affect the summation.  (If this isn't done, any sum that includes a NaN
% will also become NaN.)
%
A = isnan(matrixIn);
matrixIn(A) = 0;

%
% For each element, we have to count how many non-NaN elements went into
% the sums.  This is so we can divide by that number to get a mean.  We use
% the same matrices to do this (ie, "eL" and "eR").
%
nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;

%
% Actually taking the mean.
%
matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;
end

