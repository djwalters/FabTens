clear all;
fileID = fopen('C:\Users\David\Documents\MSU Research\Doctoral Work\Mechanical Testing\Radiation Recrystallization\PhD Work\Spheres\Images\s.batman(2).txt');

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
            StrucThickHist = [dataArray{1:end-1}];
            break
        end
    end
end
fclose(fileID);