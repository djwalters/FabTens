function [idx,PixSize,VolFrac,EigVals,EigVecs,...
    MeanStrucThick,StrucThickHist,MeanStrucSep,StrucSepHist,RootPath]...
    = CTAnData(FileIn)
% CTAnData.m
% function [idx,PixSize,VolFrac,EigVals,EigVecs,...
%     MeanStrucThick,StrucThickHist,RootPath] = CTAnData(FileIn)
%
% This function reads in the data output from CTAn for a 3D volume analysis
% with adavanced anisotropy data turned on.  This produces the necessary
% variables for visualizing the Mean Intercept Length tensor as it relates
% to material texture (not to be confused with bonding texture).
%
% INPUTS:
%   FileIn: Enter either 1 or 2 to specify whether to read in a single file
%   (1) or multiple files (2).  In the former case, this allows you to
%   directly select the file to read in.  In the latter case, this allows
%   you to select a directory containing folders marked with the time of
%   the scan which contain the appropriate .txt files. E.g.:
%       >Time Lapse:
%           >0830
%               >0830Analysis_3D.txt
%               >0830Analysis_strucThick.txt
%           >1130
%               >1130Analysis_3D.txt
%               >1130Analysis_strucThick.txt
%           >1430
%               >1430Analysis_3D.txt
%               >1430Analysis_strucThick.txt
%        etc.
%
%
% OUTPUTS:
%       idx: Elapsed time identifier used to label plots and seperate data
%       based on time during experiment.  Derived from the folder stucture.
%       See help for subfunction FileImport in this file for the
%       appropriate folder structure to use with this function.
%
%       PixSize: Outputs the pixel/voxel resolution of the CT Scan for
%       proper dimensional scaling.
%       
%       VolFrac: Ice volume fraction computed for the volume of interest.
%
%       EigVals: Eigen Values of the MIL Fabric Tensor reported from CTAn.
%
%       EigVecs: Eigen Vectors of the MIL Fabric Tensor reported from CTAn.
%
%       MeanStrucThick: The mean value of the structure thickness as
%       computed from CTAn.  This is produced by inscribing a sphere along
%       all voxels of the 3-D structures' skeleton.  The structure
%       thickness is the radius of that sphere.  In addition, a histogram
%       of the distribution of the computed structure thicknesses is
%       plotted.
%
%       StruckThickHist: Histogram data of structure thickness as computed 
%       from CTAn.
%
%       RootPath: String variable containing the directory for starting
%       file selection in the correct folder.
%
% Version: 1.0 - October 23, 2014
% AUTHOR: David J. Walters; Montana State University

%% Import Files
switch FileIn
    case 1
        LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
%         LocalPath = 'C:\Users\David\Documents\MSU Research\Doctoral Work\Mechanical Testing\Radiation Recrystallization\PhD Work\';
        [CTAnStereoFile,CTAnStereoPath] = uigetfile(...
            [LocalPath,'.txt'],'Select CTAn Stereology Data');
        CTAnFullPath = [CTAnStereoPath,CTAnStereoFile];
        
        [PixSize,VolFrac,EigVals,EigVecs] = importCTAnData(CTAnFullPath);
        EigVecs = EigVecs';
        VolFrac = VolFrac/100;
        
        LocalPath = CTAnStereoPath;
        [CTAnStrucThickFile,CTAnStrucThickPath] = uigetfile(...
        [LocalPath,'.txt'],'Select CTAn Structure Thickness Data');
        StrucThickFullPath = [CTAnStrucThickPath,CTAnStrucThickFile];
        [MeanStrucThick,StrucThickHist{1}] = ...
            importStruc(StrucThickFullPath,12);
        [CTAnStrucSepFile,CTAnStrucSepPath] = uigetfile(...
        [LocalPath,'.txt'],'Select CTAn Structure Seperation Data');
        StrucSepFullPath = [CTAnStrucSepPath,CTAnStrucSepFile];
        [MeanStrucSep,StrucSepHist{1}] = ...
            importStruc(StrucSepFullPath,12);
        idx = 0;
        endtime = idx+1;
        RootPath = LocalPath;
    case 2
        % Restrict local path for selecting data folder, and select file
        % from GUI selection
        LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
        RootPath = uigetdir(LocalPath,...
            'Select root directory of segmentation results');
        % Get times from folders of specific tests
        TimeFolds = GetFolds(RootPath)
        
        % Loop through to get path for each data set and record the time
        for i = 1:length(TimeFolds)
            FilePath = fullfile(RootPath,TimeFolds{i});
            Files = dir(fullfile(FilePath,'*.txt'));
            if size(Files,1) > 3
                error('Cannot have more than three .txt file in results time folder.  Please reduce the number of .csv files to 1 in each time folder containing results')
            end
            CTAnFullPath1 = fullfile(FilePath,Files(1).name);
            CTAnFullPath2 = fullfile(FilePath,Files(2).name);
            CTAnFullPath3 = fullfile(FilePath,Files(3).name);
            [PixSize(i),VolFrac(i),EigVals(:,i),EigVecs(:,:,i)] = importCTAnData(CTAnFullPath1);
%             EigVecs(:,:,i) = (EigVecs{i})';
            VolFrac(i) = VolFrac(i)/100;
            
            [MeanStrucThick(i),StrucThickHist{i}] = ...
                importStruc(CTAnFullPath2,12);
        
            [MeanStrucSep(i),StrucSepHist{i}] = ...
                importStruc(CTAnFullPath3,12);
        
            t1{i} = datevec(TimeFolds{i},'HHMM');
            idx(i) = etime(t1{i},t1{1})/60^2;
        end
        endtime = idx(i)+1;
        days = 1;
end

if length(idx)>1
    PixSize = PixSize(1);
end

end

function [PixSize,VolFrac,EigVals,EigVecs] = importCTAnData(filename)

%% Open the text file.
fileID = fopen(filename,'r');
delimiter = ',';

%% Import PixSize
% Initialize variables for PixSize.
startRow = 12;
endRow = 12;

% Format string for each line of text:
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%*s%*s%*s%[^\n\r]';

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end
% Create output variable
PixSize = dataArray{:, 1};
% clear dataArray
%% Import VolFrac

% Initialize variables for VolFrac.
% startRow = 17;
% endRow = 17;
startRow = 5;   %Really row 17 but using previous data 12+5=17
endRow = 5;     %Really row 17 but using previous data 12+5=17

% Format string for each line of text:
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%*s%*s%*s%[^\n\r]';

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

% Create output variable
VolFrac = [dataArray{1:end-1}];

%% Import EigVals
% Initialize variables for EigVals.
% startRow = 61;
% endRow = 61;
startRow = 44;  %Really row 28 but using previous data 17+44=61
endRow = 44;    %Really row 28 but using previous data 17+44=61

% Format string for each line of text:
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%f%f%*s%[^\n\r]';

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

% Create output variable
EigVals = [dataArray{1:end-1}];
EigVals = EigVals';

%% Import EigVecs
% Initialize variables for EigVecs.
% startRow = 91;
% endRow = 93;
startRow = 30;  %Really row 91 but using previous data 61+30=91
endRow = 32;    %Really row 93 but using previous data 61+32=93
% Format string for each line of text:
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%f%f%*s%[^\n\r]';

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Create output variable
EigVecs = [dataArray{1:end-1}];
EigVecs = EigVecs'

%% Close the text file.
fclose(fileID);
end

function [MeanStrucThick,StrucThickHist] = ...
    importStruc(filename,intervals)
%% Open the text file.
fileID = fopen(filename,'r');
delimiter = ',';

%% Import Mean Structure Thickness
% Initialize variables.
startRow = 15;
endRow = 15;

% Format string for each line of text:
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%*s%[^\n\r]';

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end
% Create output variable
MeanStrucThick = [dataArray{1:end-1}];

%% Import Historgram of Structure Thickness
% Initialize variables.
startRow = 4;
endRow = inf;

% Format string for each line of text:
%   column2: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%*s%f%[^\n\r]';

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Create output variable
StrucThickHist = [dataArray{1:end-1}];

%% Close the text file.
fclose(fileID);
end
