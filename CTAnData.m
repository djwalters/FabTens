function [idx,PixSize,VolFrac,EigVals,EigVecs,...
    MeanStrucThick,StrucThickHist] = CTAnData(FileIn)

% function [PixSize,VolFrac,EigVals,EigVecs,MeanStrucThick,StrucThickHist,StrucThickSD]...
%     = CTAnData(FileIn)
%[VF,EVal,EVec,Thick] = CTAnData(argin)
% This function reads in the data output from CTAn for a 3D volume analysis
% with adavanced anisotropy data turned on.  This produces the necessary
% variables for visualizing the Mean Intercept Length tensor as it relates
% to material texture (not to be confused with bonding texture).

[idx,PixSize,VolFrac,EigVals,EigVecs,...
    MeanStrucThick,StrucThickHist] = FileImport(FileIn);

PixSize = PixSize{1};

% LocalPath = CTAnStereoPath;
% [CTAnStrucThickFile,CTAnStrucThickPath] = uigetfile(...
%     [LocalPath,'.txt'],'Select CTAn Structure Thickness Data');
% StrucThickFullPath = [CTAnStrucThickPath,CTAnStrucThickFile];
% [MeanStrucThick,StrucThickHist,StrucThickSD] = ...
%     importStrucThick(StrucThickFullPath,12);
end

function [idx,PixSize,VolFrac,EigVals,EigVecs,...
    MeanStrucThick,StrucThickHist] = FileImport(FileIn)
% FileImport is a function which imports output from CTAn which contains
% MIL Fabric tensor values, microstructural quantities, and structure
% thickness information
%
% INPUTS:
%   FileIn: Enter either 1 or 2 to specify whether to read in a single file
%   (1) or multiple files (2).  In the former case, this allows you to
%   directly select the file to read in.  In the latter case, this allows
%   you to select a directory containing folders marked with the time of
%   the scan which contain the appropriate .txt files. E.g.:
%       >Time Lapse:
%           >0830
%               >0830Analysis.txt
%           >1130
%               >1130Analysis.txt
%           >1430
%               >1430Analysis.txt
%        etc.
%
% OUTPUTS:
%   hdr: Contains the header names of the imported data
%   bonds: Output array of the entire imported dataset
%   idx: Elapsed time of each scan.  This utilizes the name of the folder
%   to calculate time.
%   endtime: A variable that sets the maximum limit of time to be used on
%   any associated plots.

switch FileIn
    case 1
        LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
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
        [MeanStrucThick,StrucThickHist] = ...
            importStrucThick(StrucThickFullPath,12);
        idx = 0;
        endtime = idx+1;
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
            if size(Files,1) > 2
                error('Cannot have more than two .txt file in results time folder.  Please reduce the number of .csv files to 1 in each time folder containing results')
            end
            CTAnFullPath1 = fullfile(FilePath,Files(1).name);
            CTAnFullPath2 = fullfile(FilePath,Files(2).name);
            [PixSize{i},VolFrac{i},EigVals{i},EigVecs{i}] = importCTAnData(CTAnFullPath1);
            EigVecs{i} = EigVecs{i}';
            VolFrac{i} = VolFrac{i}/100;
            
            [MeanStrucThick{i},StrucThickHist{i}] = ...
            importStrucThick(CTAnFullPath2,12);
%             fid = fopen(FullPath1, 'r');  %Open read only
%             % Set headers and number of columns of data (24 columns currently)
%             hdr{i} = textscan(fid, '%s',24,'delimiter',',');
%             % Scan the numerical data (11 columns currently)
%             bonds{i} = textscan(fid...
%                 ,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'...
%                 ,'delimiter',',','HeaderLines',1);
%             fclose(fid);
            t1{i} = datevec(TimeFolds{i},'HHMM');
            idx(i) = etime(t1{i},t1{1})/60^2;
        end
        endtime = idx(i)+1;
        days = 1;
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
PixSize = [dataArray{1:end-1}];
clear dataArray
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
% startRow = 28;
% endRow = 30;
startRow = 11;  %Really row 28 but using previous data 17+11=28
endRow = 13;    %Really row 28 but using previous data 17+13=30

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
EigVals = [dataArray{1:end-1}];

%% Import EigVecs
% Initialize variables for EigVecs.
% startRow = 91;
% endRow = 93;
startRow = 61;  %Really row 91 but using previous data 30+61=91
endRow = 63;    %Really row 93 but using previous data 30+63=93
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

%% Close the text file.
fclose(fileID);
end

function [MeanStrucThick,StrucThickHist] = ...
    importStrucThick(filename,intervals)
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

% %% Import Structure Thickness Standard Deviation
% % Initialize variables.
% startRow = 1;
% endRow = 1;
% 
% % Format string for each line of text:
% %   column3: double (%f)
% % For more information, see the TEXTSCAN documentation.
% formatSpec = '%*s%*s%f%*s%[^\n\r]';
% 
% % Read columns of data according to format string.
% % This call is based on the structure of the file used to generate this
% % code. If an error occurs for a different file, try regenerating the code
% % from the Import Tool.
% textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% for block=2:length(startRow)
%     frewind(fileID);
%     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
%     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
%     dataArray{1} = [dataArray{1};dataArrayBlock{1}];
% end
% 
% % Create output variable
% StrucThickSD = [dataArray{1:end-1}];
%% Close the text file.
fclose(fileID);
end
