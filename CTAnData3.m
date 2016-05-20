function [idx,PixSize,VolFrac,EigVals,LEigVals,UEigVals,EigVecs,...
    MeanStrucThick,StrucThickHist,MeanStrucSep,StrucSepHist,RootPath]...
    = CTAnData3(FileIn)
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
% Version: 1.1 - Updated filepaths, October 29, 2015
% AUTHOR: David J. Walters; Montana State University

%% Import Files
switch FileIn
    case 1
%         LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
%         LocalPath = 'C:\Users\David\Documents\MSU Research\Doctoral Work\Mechanical Testing\Radiation Recrystallization\PhD Work\';
%         [CTAnStereoFile,CTAnStereoPath] = uigetfile(...
%             [LocalPath,'.txt'],'Select CTAn Stereology Data');
%         CTAnFullPath = [CTAnStereoPath,CTAnStereoFile];
        CTAnFullPath = 'C:\Users\David\Documents\MSU Research\Doctoral Work\Mechanical Testing\Radiation Recrystallization\PhD Work\Spheres\Images\s.batman.txt';
        [PixSize,VolFrac,EigVals,LEigVals,UEigVals,EigVecs,...
    MeanStrucThick,StrucThickHist{1},...
    MeanStrucSep,StrucSepHist{1}] = importCTAnAll(CTAnFullPath);
%         [PixSize,VolFrac,EigVals,LEigVals,UEigVals,EigVecs] = importCTAnData(CTAnFullPath);
        EigVecs = EigVecs';
        VolFrac = VolFrac/100;
        
%         LocalPath = CTAnStereoPath;
%         [CTAnStrucThickFile,CTAnStrucThickPath] = uigetfile(...
%         [LocalPath,'.txt'],'Select CTAn Structure Thickness Data');
%         StrucThickFullPath = [CTAnStrucThickPath,CTAnStrucThickFile];
%         [MeanStrucThick,StrucThickHist{1}] = ...
%             importStruc(StrucThickFullPath,12);
%         [CTAnStrucSepFile,CTAnStrucSepPath] = uigetfile(...
%         [LocalPath,'.txt'],'Select CTAn Structure Seperation Data');
%         StrucSepFullPath = [CTAnStrucSepPath,CTAnStrucSepFile];
%         [MeanStrucSep,StrucSepHist{1}] = ...
%             importStruc(StrucSepFullPath,12);
        idx = 0;
        endtime = idx+1;
%         RootPath = LocalPath;
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
%             CTAnFullPath2 = fullfile(FilePath,Files(2).name);
%             CTAnFullPath3 = fullfile(FilePath,Files(3).name);
            [PixSize(i),VolFrac(i),EigVals(:,i),...
                LEigVals(:,i),UEigVals(:,i),EigVecs(:,:,i),...
    MeanStrucThick(i),StrucThickHist{i},...
    MeanStrucSep(i),StrucSepHist{i}] = importCTAnAll(CTAnFullPath1);
%             [PixSize(i),VolFrac(i),EigVals(:,i),LEigVals(:,i),UEigVals(:,i),EigVecs(:,:,i)] = importCTAnData(CTAnFullPath1);
%             EigVecs(:,:,i) = (EigVecs{i})';
            VolFrac(i) = VolFrac(i)/100;
            
%             [MeanStrucThick(i),StrucThickHist{i}] = ...
%                 importStruc(CTAnFullPath2,12);
%         
%             [MeanStrucSep(i),StrucSepHist{i}] = ...
%                 importStruc(CTAnFullPath3,12);
        
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
function [PixSize,VolFrac,EigVals,LEigVals,UEigVals,EigVecs,...
    MeanStrucThick,StrucThickHist,...
    MeanStrucSep,StrucSepHist] = importCTAnAll(filename)
fileID = fopen(filename);

% Open function
counter = 0;                                % line counter (sloppy but works)

%% Get Pixel Size
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Pixel size,,'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            %             break                           % we found it, lets go home
            formatSpec = '%*s%*s%f%*s%[^\n\r]';
            delimiter = ',';
            % Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            %             textscan(fileID, '%[^\n\r]', 'ReturnOnError', false)
            dataArray = textscan(tline, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            PixSize = dataArray{:, 1};
            break
        end
    end
end
%% Get Volume Fraction
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Percent object volume,'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            %             break                           % we found it, lets go home
            formatSpec = '%*s%*s%f%*s%[^\n\r]';
            delimiter = ',';
            % Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            %             textscan(fileID, '%[^\n\r]', 'ReturnOnError', false)
            dataArray = textscan(tline, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            VolFrac = dataArray{:,1};
            break
        end
    end
end
%% Get Eigen Values
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Principal MILs,'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            %             break                           % we found it, lets go home
            formatSpec = '%*s%*s%f%f%f%*s%[^\n\r]';
            delimiter = ',';
            % Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            %             textscan(fileID, '%[^\n\r]', 'ReturnOnError', false)
            dataArray = textscan(tline, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            EigVals = [dataArray{1:end-1}];
            EigVals = EigVals';
            break
        end
    end
end
%% Get Eigen Value Confidence Intervals
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, '95% Confidence intervals for principals:'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            %             break                           % we found it, lets go home
            for i = 1:3
                tline = fgetl(fileID);
                counter = counter + 1;
                
                formatSpec = '%7f%*9c%7f%[^\n\r]';
                %             delimiter = ',';
                % Read columns of data according to format string.
                % This call is based on the structure of the file used to generate this
                % code. If an error occurs for a different file, try regenerating the code
                % from the Import Tool.
                dataArray = textscan(tline, formatSpec);
                [LEigVals(i)] = dataArray{1, 1};
                [UEigVals(i)] = dataArray{1, 2};
            end
            
            break
        end
    end
end
%% Get Eigen Value Confidence Intervals
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'E-vector 1,'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            %             break                           % we found it, lets go home
            for i = 1:3
                formatSpec = '%*s%*s%f%f%f%*s%[^\n\r]';
                delimiter = ',';
                % Read columns of data according to format string.
                % This call is based on the structure of the file used to generate this
                % code. If an error occurs for a different file, try regenerating the code
                % from the Import Tool.
                %             textscan(fileID, '%[^\n\r]', 'ReturnOnError', false)
                dataArray = textscan(tline, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
                [EigVecs(i,:)] = [dataArray{1:end-1}]
                tline = fgetl(fileID);
                counter = counter + 1;
                
            end
            EigVecs = EigVecs'
            
            break
        end
    end
end
%% Get Mean Structure Thickness
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Structure thickness,St.Th,'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            %             break                           % we found it, lets go home
            formatSpec = '%*s%*s%f%*s%[^\n\r]';
            delimiter = ',';
            % Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            %             textscan(fileID, '%[^\n\r]', 'ReturnOnError', false)
            dataArray = textscan(tline, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            MeanStrucThick = [dataArray{1:end-1}];
            break
        end
    end
end
%% Get Structure Thickness Histogram
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Structure thickness distribution,St.Th'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            startRow = counter+3;
            continue
        end
    end
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Standard deviation of structure thickness'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            endRow = counter-1;
            
            % Format string for each line of text:
            %   column2: double (%f)
            %	column4: double (%f)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%*s%f%*s%f%[^\n\r]';
            delimiter = ',';
            % Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            frewind(fileID)
            textscan(fileID, '%[^\n\r]', startRow(1),'Delimiter','\n', 'ReturnOnError', false);
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1), 'Delimiter', delimiter, 'ReturnOnError', false);
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
            break
        end
    end
end
counter = endRow - 1;
%% Get Mean Structure Seperation
while 1                                     % infinite loop
    tline = fgetl(fileID)                     % read a line
    counter = counter + 1                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Structure separation,St.Sp,'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            %             break                           % we found it, lets go home
            formatSpec = '%*s%*s%f%*s%[^\n\r]';
            delimiter = ',';
            % Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            %             textscan(fileID, '%[^\n\r]', 'ReturnOnError', false)
            dataArray = textscan(tline, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            MeanStrucSep = [dataArray{1:end-1}];
            break
        end
    end
end
%% Get Structure Seperation Histogram
while 1                                     % infinite loop
    tline = fgetl(fileID);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Structure separation distribution,St.Sp'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            startRow = counter+3
            continue
        end
    end
    if ischar(tline)                        % if the line is string
        U = strfind(tline, 'Standard deviation of structure separation'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            endRow = counter-1
            
            % Format string for each line of text:
            %   column2: double (%f)
            %	column4: double (%f)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%*s%f%*s%f%[^\n\r]';
            delimiter = ',';
            % Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            frewind(fileID)
            textscan(fileID, '%[^\n\r]', startRow(1),'Delimiter','\n', 'ReturnOnError', false)
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1), 'Delimiter', delimiter, 'ReturnOnError', false);
            for block=2:length(startRow)
                frewind(fileID);
                textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end
            
            % Create output variable
            StrucSepHist = [dataArray{1:end-1}];
            break
        end
    end
end
fclose(fileID);
end

% % function [PixSize,VolFrac,EigVals,LEigVals,UEigVals,EigVecs] = importCTAnData(filename)
% % 
% % %% Open the text file.
% % fileID = fopen(filename,'r');
% % delimiter = ',';
% % 
% % %% Import PixSize
% % % Initialize variables for PixSize.
% % startRow = 12;
% % endRow = 12;
% % 
% % % Format string for each line of text:
% % %   column3: double (%f)
% % % For more information, see the TEXTSCAN documentation.
% % formatSpec = '%*s%*s%f%*s%*s%*s%[^\n\r]';
% % 
% % % Read columns of data according to format string.
% % % This call is based on the structure of the file used to generate this
% % % code. If an error occurs for a different file, try regenerating the code
% % % from the Import Tool.
% % textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
% % dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% % for block=2:length(startRow)
% %     frewind(fileID);
% %     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
% %     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% %     dataArray{1} = [dataArray{1};dataArrayBlock{1}];
% % end
% % % Create output variable
% % PixSize = dataArray{:, 1};
% % % clear dataArray
% % %% Import VolFrac
% % 
% % % Initialize variables for VolFrac.
% % % startRow = 17;
% % % endRow = 17;
% % startRow = 5;   %Really row 17 but using previous data 12+5=17
% % endRow = 5;     %Really row 17 but using previous data 12+5=17
% % 
% % % Format string for each line of text:
% % %   column3: double (%f)
% % % For more information, see the TEXTSCAN documentation.
% % formatSpec = '%*s%*s%f%*s%*s%*s%[^\n\r]';
% % 
% % % Read columns of data according to format string.
% % % This call is based on the structure of the file used to generate this
% % % code. If an error occurs for a different file, try regenerating the code
% % % from the Import Tool.
% % textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
% % dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% % for block=2:length(startRow)
% %     frewind(fileID);
% %     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
% %     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% %     dataArray{1} = [dataArray{1};dataArrayBlock{1}];
% % end
% % 
% % % Create output variable
% % VolFrac = [dataArray{1:end-1}];
% % 
% % %% Import EigVals
% % % Initialize variables for EigVals.
% % % startRow = 61;
% % % endRow = 61;
% % startRow = 44;  %Really row 28 but using previous data 17+44=61
% % endRow = 44;    %Really row 28 but using previous data 17+44=61
% % % startRow = 40;
% % % endRow = 40;
% % % Format string for each line of text:
% % %   column3: double (%f)
% % %	column4: double (%f)
% % %   column5: double (%f)
% % % For more information, see the TEXTSCAN documentation.
% % formatSpec = '%*s%*s%f%f%f%*s%[^\n\r]';
% % 
% % % Read columns of data according to format string.
% % % This call is based on the structure of the file used to generate this
% % % code. If an error occurs for a different file, try regenerating the code
% % % from the Import Tool.
% % textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
% % dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% % for block=2:length(startRow)
% %     frewind(fileID);
% %     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
% %     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% %     dataArray{1} = [dataArray{1};dataArrayBlock{1}];
% % end
% % 
% % % Create output variable
% % EigVals = [dataArray{1:end-1}];
% % EigVals = EigVals';
% % 
% % %% Import Bounds on EigenVals
% % % Initialize variables.
% %     startRow = 13;
% %     endRow = 15;
% % 
% % 
% % %% Format string for each line of text:
% % %   column2: text (%s)
% % %	column4: double (%f)
% % % For more information, see the TEXTSCAN documentation.
% % formatSpec = '%7f%11*s%7f%[^\n\r]';
% % 
% % %% Read columns of data according to format string.
% % % This call is based on the structure of the file used to generate this
% % % code. If an error occurs for a different file, try regenerating the code
% % % from the Import Tool.
% % textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
% % dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
% % for block=2:length(startRow)
% %     frewind(fileID);
% %     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
% %     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
% %     for col=1:length(dataArray)
% %         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
% %     end
% % end
% % 
% % 
% % %% Post processing for unimportable data.
% % % No unimportable data rules were applied during the import, so no post
% % % processing code is included. To generate code which works for
% % % unimportable data, select unimportable cells in a file and regenerate the
% % % script.
% % 
% % %% Allocate imported array to column variable names
% % LEigVals = dataArray{:, 1};
% % UEigVals = dataArray{:, 2};
% % 
% % 
% % %% Import EigVecs
% % % Initialize variables for EigVecs.
% % % startRow = 91;
% % % endRow = 93;
% % startRow = 15;  %Really row 91 but using previous data 61+30=91
% % endRow = 17;    %Really row 93 but using previous data 61+32=93
% % % startRow = 29;
% % % endRow = 31;
% % % Format string for each line of text:
% % %   column3: double (%f)
% % %	column4: double (%f)
% % %   column5: double (%f)
% % % For more information, see the TEXTSCAN documentation.
% % formatSpec = '%*s%*s%f%f%f%*s%[^\n\r]';
% % 
% % % Read columns of data according to format string.
% % % This call is based on the structure of the file used to generate this
% % % code. If an error occurs for a different file, try regenerating the code
% % % from the Import Tool.
% % textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
% % dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% % for block=2:length(startRow)
% %     frewind(fileID);
% %     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
% %     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
% %     for col=1:length(dataArray)
% %         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
% %     end
% % end
% % 
% % % Create output variable
% % EigVecs = [dataArray{1:end-1}];
% % EigVecs = EigVecs';
% % 
% % %% Close the text file.
% % fclose(fileID);
% % end
% 
% function [MeanStrucThick,StrucThickHist] = ...
%     importStruc(filename,intervals)
% %% Open the text file.
% fileID = fopen(filename,'r');
% delimiter = ',';
% 
% %% Import Mean Structure Thickness
% % Initialize variables.
% startRow = 15;
% endRow = 15;
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
% % Create output variable
% MeanStrucThick = [dataArray{1:end-1}];
% 
% %% Import Historgram of Structure Thickness
% % Initialize variables.
% startRow = 4;
% endRow = inf;
% 
% % Format string for each line of text:
% %   column2: double (%f)
% %	column4: double (%f)
% % For more information, see the TEXTSCAN documentation.
% formatSpec = '%*s%f%*s%f%[^\n\r]';
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
%     for col=1:length(dataArray)
%         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%     end
% end
% 
% % Create output variable
% StrucThickHist = [dataArray{1:end-1}];
% 
% %% Close the text file.
% fclose(fileID);
% end
