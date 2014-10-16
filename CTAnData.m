function [PixSize,VolFrac,EigVals,EigVecs,MeanStrucThick,StrucThickHist,StrucThickSD]...
    = CTAnData(~)
%[VF,EVal,EVec,Thick] = CTAnData(argin)
% This function reads in the data output from CTAn for a 3D volume analysis
% with adavanced anisotropy data turned on.  This produces the necessary
% variables for visualizing the Mean Intercept Length tensor as it relates
% to material texture (not to be confused with bonding texture).

LocalPath = 'C:\Doctoral Researach\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Linux FTP Mirror\PhD Work\MicroMechanics\Matlab 3D Segmentation Results\';
[CTAnStereoFile,CTAnStereoPath] = uigetfile(...
    [LocalPath,'.txt'],'Select CTAn Stereology Data');
CTAnFullPath = [CTAnStereoPath,CTAnStereoFile];

[PixSize,VolFrac,EigVals,EigVecs] = importCTAnData1(CTAnFullPath);
EigVecs = EigVecs';
VolFrac = VolFrac/100;

LocalPath = CTAnStereoPath;
[CTAnStrucThickFile,CTAnStrucThickPath] = uigetfile(...
    [LocalPath,'.txt'],'Select CTAn Structure Thickness Data');
StrucThickFullPath = [CTAnStrucThickPath,CTAnStrucThickFile];
[MeanStrucThick,StrucThickHist,StrucThickSD] = ...
    importStrucThick(StrucThickFullPath,12);
end

function [PixSize,VolFrac,EigVals,EigVecs] = importCTAnData1(filename)
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

function [MeanStrucThick,StrucThickHist,StrucThickSD] = ...
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
endRow = startRow + (intervals-1);

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

%% Import Structure Thickness Standard Deviation
% Initialize variables.
startRow = 1;
endRow = 1;

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
StrucThickSD = [dataArray{1:end-1}];

%% Close the text file.
fclose(fileID);
end
