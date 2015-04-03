function [F2,F2Ci,F4,F4Ci,bondRad,grainRad,meanBondRad,meanGrainRad,...
    coordNum,shapeFac,idx,endtime,spatialLabel]...
    = ContactTensor(FileIn,BondVec,Plot,varargin)
% ContactTensor calculates the coefficients for the 2nd and 4th order
% contact tensor derived from 3-D segmentation data.  For 2-D stereology,
% refer to Dr. Rich Shertzer's code for analysis.
% Code also calculates confidence intervals via a bootstrap sampling method
% and also analyzes the statistical significance of the 2nd order tensor
% vs. the 4th order tensor.
%
%   INPUTS: 
%       FileIn: specifies whether a single file is used or multiple files
%       chronologically ordered. 1 for single file, 2 for multiple files
%
%       BondVec: String input specifying method for calculating contact
%       tensor; options
%           'planenorm' - uses bond plane normal vector directly.  Caution
%           must be applied when using this value as bonds must be
%           sufficiently large to accurately define the plane.
%
%           'geocenter2center' - utilizes the vector connecting each grain
%           center (based on center of effective spherical grain radius) 
%           in the grain-bond pair.  Can be used when 'planenorm' is not 
%           giving accurate results
%
%           'masscenter2center' - utilizes the vector connecting each grain
%           center (mass centroid) in the grain-bond pair.  Can be used 
%           when 'planenorm' is not giving accurate results.
%
%       Plot: Indicate 1 to include plots, 0 to exclude plots
%
%       varargin{1}: Optional input to specify the spatial pixel dimension in
%       micrometers (um).  If left blank, all spatial dimensions will be
%       based on per pixel rather than per micrometer.
%
%       varargin{2}: Optional input of the rootpath directory for starting
%       file selection in the correct folder.  This is used when function
%       is called in a higher level function where file selection is
%       already taking place.
%
%   OUTPUS: 
%       F2(i,j,n): 2nd order contact tensor.  i and j reference the matrix
%       position within an individual tensor (i=1,2,3;j=1,2,3) and n
%       specifies the volume in chronological order with which the contact
%       tensor is associated.
%       
%       F2Ci(i,j,n,b): Upper and lower bonds of the 95% confidence interval of
%       2nd order contact tensor.  i, j, and n are the same as above.  b 
%       indicates the upper or lower bounds of the confidence interval. 1
%       indicates the lower bound, 2 indicates the upper bound.
%
%       F4(i,j,k,l,n): 4th order contact tensor. i,j,k,l are the indices of
%       the 4th order contact tensor (i,j,k,l=1,2,3). n specifies the 
%       volume in chronological order with which the contact tensor is
%       associated.
%
%       F4Ci(i,j,k,l,n,b): Upper and lower bonds of the 95% confidence 
%       interval of 2nd order contact tensor.  i,j,k,l,n are the same as 
%       above.  b indicates the upper or lower bounds of the confidence 
%       interval. 1 indicates the lower bound, 2 indicates the upper bound.
%
%       bondRad{n}: A vector of the distribution of bond radii contained
%       within a cell structure for time series data.
%
%       grainRad{n}: Same as bondRad except the distribution of grain
%       radii.
%
%       meanBondRad: Mean 3-D bond radius.  This assumes a circular contact
%       since it is based off the calculated bond area.
%
%       meanGrainRad: Mean 3-D grain radius based on the largest inscribed
%       sphere of each individual grain identified in 3-D segmentation
%
%       coordNum: Mean 3-D coordination number.  This number represents the
%       average number of contacts per grain.
%
%       shapeFac: A calculation of the scalar shape factor.  This is the
%       ratio of the spherical grain gradius (meanGrainRad) to the distance
%       of the grain center to bond center.
%
%       idx: Elapsed time identifier used to label plots and seperate data
%       based on time during experiment.  Derived from the folder stucture.
%       See help for subfunction FileImport in this file for the
%       appropriate folder structure to use with this function.
%
%       endtime: Variable listing the final elapsed time used in plotting
%       functions
%
%       spatialLabel: String variable set here to use in plotting functions
%       in other scripts and functions.
%
%   EXTERNAL FUNCTION CALLS:
%       MnDev.m
%       ConvgTest.m
%       padcat.m
%       errorbarwidth.m
%       GetFolds.m
%       
%
% Author: David J. Walters, Montana State University
%
% Version 1.0 - Calculates contact tensor for a single file from 3-D
% segmentation data.  May 2, 2014
%
% Version 1.1 - Allows the import of timelapsed data from a directory of
% multiple segmentation data sets.  May 15, 2014
%
% Version 1.1.1 - Plots grain radius and bond radius distribution.
%
% Version 1.2 - Includes calculation of grain shape factor from the
% difference between spherical grain radius and the distance from spherical
% grain center to the bond center.  June 17, 2014
%
% Version 1.2.5 - Removed plots to seperate external function called
% ContactPlot.m.  October 22, 2014
%
% Version 1.3.0 - Added switch to how bond orientations are calculated.
% April 2, 2015

%% Check Inputs
if isempty(varargin)
    pixdim = 1;
    spatialLabel = '( pix )';
    fprintf('All spatial dimensions are in pixels since no pixel dimension was provided.\n')
    RootPath = 0;
elseif length(varargin) == 1
    pixdim = varargin{1};
    spatialLabel = '( \mum )';
    RootPath = 0;
elseif length(varargin) == 2
    pixdim = varargin{1};
    spatialLabel = '( \mum )';
    RootPath = varargin{2};
else
    error('Too many input arguments')
end
%% Initialize file input cases for either single analysis or a series of analyses
nboot = 1000;
[hdr,bonds,idx,endtime] = FileImport(FileIn,RootPath);   %Subfunction for importing files

% Create a waitbar to show progress of computation including a cancel
% button
h = waitbar(0,'Please wait...',...
    'Name','Calculating Contact Tensor...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)
waitprev = 0;
for n = 1:length(idx)
    sprintf('idx(n)= %2.1f',idx(n))
    % Build data structure of output
    for i = 1:length(hdr{n}{1});
        % Searches for numerical position of the end of character string
        % before unit definition, i.e. rb(um), idx = 3 results in (data.rb)
        label = strfind(hdr{n}{1}{i},'(');
        data{n}.(hdr{n}{1}{i}(1:label-1)) = bonds{n}{i};
    end
    
    switch BondVec
        case 'planenorm'
            % Generates bond plane cell array directly from bond plane
            % normal vectors calculated during 3D segmentation
            BondPlanes{n} = [data{n}.BNx,data{n}.BNy,data{n}.BNz];
            
        case 'geocenter2center'
            % Generates bond plane cell array from the coordinates of the
            % centers of each grain in a grain-bond pair based on the
            % geometric center (grain radius center)
            
            % Calculate center to center vectors for all data at once
            BP = [data{n}.G2RCx - data{n}.G1RCx,...
                data{n}.G2RCy - data{n}.G1RCy,...
                data{n}.G2RCz - data{n}.G1RCz];
            
            %Need unit vectors.  These are generated here, one by one
            for i = 1:size(BP,1)
                BondPlanes{n}(i,:) = BP(i,:)/norm(BP(i,:));
            end
        case 'masscenter2center'
            % Generates bond plane cell array from the coordinates of the
            % centers of each grain in a grain-bond pair based on the
            % mass center (mass centroid)
            BP = [data{n}.G2Cx - data{n}.G1Cx,...
                data{n}.G2Cy - data{n}.G1Cy,...
                data{n}.G2Cz - data{n}.G1Cz];
            for i = 1:size(BP,1)
                BondPlanes{n}(i,:) = BP(i,:)/norm(BP(i,:));
            end
    end
            
    
    % Bootstrap sampling for confidence intervals
    for p = 1:nboot
        % Check for Cancel button press
        if getappdata(h,'canceling')
            delete(h)
            error('Canceled')
        end
        waitval = waitprev+((p/nboot)/3)/length(idx);
        waitbar(waitval,h,...
            sprintf('Bootstrap sampling in progress... Total progress %12.1f%%',waitval*100))
        % Resample directional data before calculating tensor coefficients
        % Creates a random sample of bond numbers that is the has the
        % same length as the number of bonds, and selected with
        % replacement.
        SampRows = randsample(data{n}.Row,length(data{n}.Row),true);
        
        % Creates array of random sampled bond orientations based on
        % bonds identified in SampRows.
        for i = 1:length(SampRows)
            BPSamp(i,:) = BondPlanes{n}(SampRows(i),:);
        end
        
        F2Samp = zeros(3,3);
        F4Samp = zeros(3,3,3,3);
        for i = 1:3
            for j = 1:3
                for m=1:size(BPSamp,1)
                    F2Samp(i,j) = F2Samp(i,j) + BPSamp(m,i)*BPSamp(m,j);
                end
                for k = 1:3
                    for l = 1:3
                        for m = 1:size(BPSamp,1)
                            F4Samp(i,j,k,l) = F4Samp(i,j,k,l) + BPSamp(m,i)...
                                *BPSamp(m,j)*BPSamp(m,k)...
                                *BPSamp(m,l);
                        end
                    end
                end
            end
        end
        F2Samp = 1/length(data{n}.Bond)*F2Samp;
        F4Samp = 1/length(data{n}.Bond)*F4Samp;
        % Store each iteration of tensor coefficients for confidence
        % intervals
        F2TotSamp(:,:,p) = F2Samp;
        F4TotSamp(:,:,:,:,p) = F4Samp;

    end
    
    %Calculate the 95% confidence interval of the resampled tensor
    %coefficients.
    F2Ci(:,:,n,:) = prctile(F2TotSamp,[5 95],3);
    F4Ci(:,:,:,:,n,:) = prctile(F4TotSamp,[5 95],5);
    clear F2Samp
    clear F4Samp
    clear BPSamp
    clear SampRows
    clear F2TotSamp 
    clear F4TotSamp
    % Calculate the 2nd order and 4th order contact tensors for the
    % entire sample of bonds.
    F2(:,:,n) = zeros(3,3);
    F4(:,:,:,:,n) = zeros(3,3,3,3);
    waitval2 = waitval;
    for i = 1:3
        waitval2 = waitval;
        for j = 1:3
            if getappdata(h,'canceling')
                delete(h)
                error('Canceled')
            end
            waitval = waitval2+(2/3)*(j/(3^2))/length(idx);
            waitbar(waitval,h,...
                sprintf('Calculating contact tensors... Total progress %12.1f%%',waitval*100))
            for m=1:size(BondPlanes{n},1)
                F2(i,j,n) = F2(i,j,n) + BondPlanes{n}(m,i)*BondPlanes{n}(m,j);
            end
            for k = 1:3
                for l = 1:3
                    for m = 1:size(BondPlanes{n},1)
                        F4(i,j,k,l,n) = F4(i,j,k,l,n) + BondPlanes{n}(m,i)...
                            *BondPlanes{n}(m,j)*BondPlanes{n}(m,k)...
                            *BondPlanes{n}(m,l);
                    end
                end
            end
        end
    end
    
    F2(:,:,n) = 1/length(data{n}.Bond)*F2(:,:,n);
    F4(:,:,:,:,n) = 1/length(data{n}.Bond)*F4(:,:,:,:,n);
    if getappdata(h,'canceling')
        delete(h)
        error('Canceled')
    end
    waitval = waitprev+(7/8)/length(idx);
    waitbar(waitval,h,...
        sprintf('Finished calculating contact tensors... Total progress %12.1f%%',waitval*100))
    
    %% Calculate deviatoric part of contact tensors - required for N2 vs N4
    % Check for Cancel button press
    [D2,D4] = MnDev(F2(:,:,n),F4(:,:,:,:,n));   % External function
    
    %% Convergence of series approximations - statistical check
    alpha = 0.05;   % Significance of test - "95%"
    [zc,z2,z4] = ConvgTest(length(data{n}.Bond),3,D2,D4,alpha);   % External Function
    % Check for Cancel button press
    if getappdata(h,'canceling')
        delete(h)
        error('Canceled')
    end
    waitval = waitprev+(1/length(idx));
    waitbar(waitval,h,...
        sprintf('Testing significance... Total progress %12.1f%%',waitval*100))
    waitprev = waitval;
    
    %% Calculate mean coordination number
    numBonds = length(data{n}.Bond);
    numGrains = length(unique([data{n}.G1 data{n}.G2]));
    coordNum(n) = 2*numBonds/numGrains;
    
    %% Calculate mean bond radius to grain radius ratio
    bondRad{n} = (data{n}.BA./pi).^0.5;
    grainRad{n} = data{n}.G1Rad;
    
    meanBondRad(n) = mean(bondRad{n});
    meanGrainRad(n) = mean(grainRad{n});
    
    %% Calculate non-spherical grain shape factor
    % Calculate the vector from the bond center to the grain radius center
    % of grain 1 of the bond.
    shapeVec{n} = [data{n}.BCx - data{n}.G1RCx,...
        data{n}.BCy - data{n}.G1RCy,...
        data{n}.BCz - data{n}.G1RCz];
    
    % Norm calculates the distance/magnitude of the vector.
    for i = 1:size(shapeVec{n},1)
        bondDistInd{n}(i) = norm(shapeVec{n}(i,:),2);
    end
    
    % Calculate the shape factor
    meanBondDist(n) = mean(bondDistInd{n});
    shapeFac(n) = meanGrainRad(n)/meanBondDist(n);
    

end

% Numbers of identified bonds and grains may vary if comparing more than
% one test.  Padcat is a function that pads arrays with "NaNs" when
% concatenating arrays together of different sizes.  The logic statements
% that follow allow for up to six analyses to be compared at once.  If more
% are needed, add more elseif statements following the syntax pattern found
% below.
if n == 1
    bondRad = bondRad{n};
    grainRad = grainRad{n};
elseif n == 2
    bondRad = padcat(bondRad{n-1},bondRad{n});
    grainRad = padcat(grainRad{n-1},grainRad{n});
elseif n == 3
    bondRad = padcat(bondRad{n-2},bondRad{n-1},bondRad{n});
    grainRad = padcat(grainRad{n-2},grainRad{n-1},grainRad{n});
elseif n == 4
    bondRad = padcat(bondRad{n-3},bondRad{n-2},bondRad{n-1},bondRad{n});
    grainRad = padcat(grainRad{n-3},grainRad{n-2},grainRad{n-1},grainRad{n});    
elseif n == 5
    bondRad = padcat(bondRad{n-4},bondRad{n-3},bondRad{n-2},bondRad{n-1},bondRad{n});
    grainRad = padcat(grainRad{n-4},grainRad{n-3},grainRad{n-2},grainRad{n-1},grainRad{n});
elseif n == 6
    bondRad = padcat(bondRad{n-5},bondRad{n-4},bondRad{n-3},bondRad{n-2},bondRad{n-1},bondRad{n});
    grainRad = padcat(grainRad{n-5},grainRad{n-4},grainRad{n-3},grainRad{n-2},grainRad{n-1},grainRad{n});    
end

% Deletes the waitbar when complete
delete(h)

%% Plots & Figures

if Plot == 1
ContactPlot(F2,F2Ci,bondRad,grainRad,...
    meanBondRad,meanGrainRad,coordNum,shapeFac,...
    idx,pixdim,endtime,spatialLabel);
end
% % Plot tensor coefficients
% ebwidth = 0.1;  % Errorbar cap width
% font = 'Palatino Linotype';
% fsize = 11;
% msize = 5;
% figure('Name','3-D Tensor Diagonals','NumberTitle','off')
% hE = errorbar(idx,F2(1,1,:),F2(1,1,:)-F2Ci(1,1,:,1),...
%     F2Ci(1,1,:,2)-F2(1,1,:),'rsquare','MarkerSize',msize);
% errorbarwidth(hE,ebwidth);
% hold on
% hE = errorbar(idx,F2(2,2,:),F2(2,2,:)-F2Ci(2,2,:,1),...
%     F2Ci(2,2,:,2)-F2(2,2,:),'bdiamond','MarkerSize',msize+1);
% errorbarwidth(hE,ebwidth);
% hold on
% hE = errorbar(idx,F2(3,3,:),F2(3,3,:)-F2Ci(3,3,:,1),...
%     F2Ci(3,3,:,2)-F2(3,3,:),'go','MarkerSize',msize);
% errorbarwidth(hE,ebwidth);
% grid
% axis([-1 (endtime) 0 1])
% y1 = ylabel('3-D Tensor Diagonals (-)');
% x1 = xlabel('Elapsed Time ( hrs )');
% set([y1 x1],'FontName',font,'FontSize',fsize)
% legend('\it{x}\rm{_1}','\it{x}\rm{_2}','\it{x}\rm{_3}',...
%     'Location','Northwest','Orientation','horizontal');
% set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)
% 
% % Plot distribution of bond radii
% figure('Name','Bond Radii Distribution','NumberTitle','off')
% hist(bondRad*pixdim,20)
% y1 = ylabel('Number of Bonds');
% x1 = xlabel(['\rho: 3-D Bond Radius ',spatialLabel]);
% set([y1 x1],'FontName',font,'FontSize',fsize)
% 
% % Plot distribution of grain radii
% figure('Name','Grain Radii Distribution','NumberTitle','off')
% hist(grainRad*pixdim,20)
% y1 = ylabel('Number of Grains');
% x1 = xlabel(['R: 3-D Grain Radius ',spatialLabel]);
% set([y1 x1],'FontName',font,'FontSize',fsize)
% 
% %Mean Bond and Grain Radii
% figure('Name','Mean Bond and Grain Radii','NumberTitle','off')
% ax(1) = subplot(1,3,1);
% plot(idx,meanBondRad*pixdim,'ko','MarkerSize',msize);
% grid
% axis([-1 (endtime) 0 200])
% y1(1) = ylabel(['\rho: Mean 3-D Bond Radius ',spatialLabel]);
% 
% ax(2) = subplot(1,3,2);
% plot(idx,meanGrainRad*pixdim,'ko','MarkerSize',msize);
% grid
% axis([-1 (endtime) 0 200])
% y1(2) = ylabel(['R: Mean 3-D Grain Radius ',spatialLabel]);
% x1 = xlabel('Elapsed Time (hrs)');
% 
% ax(3) = subplot(1,3,3);
% plot(idx,meanBondRad./meanGrainRad,'ko','MarkerSize',msize);
% grid
% axis([-1 (endtime) 0 1])
% y1(3) = ylabel('Mean 3-D Bond-to-Grain Radius Ratio (-)');
% set([y1 x1],'FontName',font,'FontSize',fsize)
% set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)
% 
% %3-D Coordination Number
% figure('Name','Mean 3-D Coordination Number','NumberTitle','off')
% plot(idx,coordNum,'ko','MarkerSize',msize);
% grid
% y1 = ylabel('Mean 3-D Coordination Number (-)');
% x1 = xlabel('Elapsed Time (hrs)');
% set([y1 x1],'FontName',font,'FontSize',fsize)
% % axis([-1 (endtime) 1 3.5])
% xlim([-1 (endtime)])
% set(gca,'XTick',0:3:endtime)
end

function [hdr,bonds,idx,endtime] = FileImport(FileIn,RootPath)
% FileImport is a function which imports 3-D segmentation data as produced
% by the Matlab function Segmentation3D run in a Linux environment.  This
% produces a .csv file containing bond and grain size, position, and
% orientation information.
%
% INPUTS:
%   FileIn: Enter either 1 or 2 to specify whether to read in a single file
%   (1) or multiple files (2).  In the former case, this allows you to
%   directly select the file to read in.  In the latter case, this allows
%   you to select a directory containing folders marked with the time of
%   the scan which contain the appropriate .csv files. E.g.:
%       >Time Lapse:
%           >0830
%               >0830Analysis.csv
%           >1130
%               >1130Analysis.csv
%           >1430
%               >1430Analysis.csv
%        etc.
%
% OUTPUTS:
%   hdr: Contains the header names of the imported data
%   bonds: Output array of the entire imported dataset
%   idx: Elapsed time of each scan.  This utilizes the name of the folder
%   to calculate time.
%   endtime: A variable that sets the maximum limit of time to be used on
%   any associated plots.
if RootPath == 0
    switch FileIn
        case 1
            LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
            [FileName,FilePath] = uigetfile([LocalPath,'.csv'],'Select Segmentation Data');
            FullPath = [FilePath,FileName];
            % Read files of segmentation data
            fid = fopen(FullPath, 'r');  %Open read only
            % Set headers and number of columns of data (24 columns currently)
            hdr{1} = textscan(fid, '%s',25,'delimiter',',');
            % Scan the numerical data (11 columns currently)
            bonds{1} = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'...
                ,'delimiter',',','HeaderLines',1);
            fclose(fid);
            idx = 0;
            endtime = idx+1;
        case 2
            % Restrict local path for selecting data folder, and select file
            % from GUI selection
            LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
            RootPath = uigetdir(LocalPath,...
                'Select root directory of segmentation results');
            % Get times from folders of specific tests
            TimeFolds = GetFolds(RootPath);
            
            % Loop through to get path for each data set and record the time
            for i = 1:length(TimeFolds)
                FilePath = fullfile(RootPath,TimeFolds{i});
                Files = dir(fullfile(FilePath,'*.csv'));
                if size(Files,1) > 1
                    error('Cannot have more than one .csv file in results time folder.  Please reduce the number of .csv files to 1 in each time folder containing results')
                end
                FullPath = fullfile(FilePath,Files.name);
                fid = fopen(FullPath, 'r');  %Open read only
                % Set headers and number of columns of data (24 columns currently)
                hdr{i} = textscan(fid, '%s',25,'delimiter',',');
                % Scan the numerical data (11 columns currently)
                bonds{i} = textscan(fid...
                    ,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'...
                    ,'delimiter',',','HeaderLines',1);
                fclose(fid);
                t1{i} = datevec(TimeFolds{i},'HHMM');
                idx(i) = etime(t1{i},t1{1})/60^2;
            end
            endtime = idx(i)+1;
            days = 1;
    end
else
    switch FileIn
        case 1
%             LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
            [FileName,FilePath] = uigetfile([RootPath,'.csv'],'Select Segmentation Data');
            FullPath = [FilePath,FileName];
%             FullPath = [RootPath,'Spheres.csv'];
            % Read files of segmentation data
            fid = fopen(FullPath, 'r');  %Open read only
            % Set headers and number of columns of data (24 columns currently)
            hdr{1} = textscan(fid, '%s',25,'delimiter',',');
            % Scan the numerical data (11 columns currently)
            bonds{1} = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'...
                ,'delimiter',',','HeaderLines',1);
            fclose(fid);
            idx = 0;
            endtime = idx+1;
        case 2
            % Restrict local path for selecting data folder, and select file
            % from GUI selection
%             LocalPath = 'C:\Doctoral Research\Mechanical Testing\Radiation Recrystallization\Fabric Tensor and ANSYS\Matlab 3D Segmentation Results\';
%             RootPath = uigetdir(LocalPath,...
%                 'Select root directory of segmentation results');
            % Get times from folders of specific tests
            TimeFolds = GetFolds(RootPath)
            
            % Loop through to get path for each data set and record the time
            for i = 1:length(TimeFolds)
                FilePath = fullfile(RootPath,TimeFolds{i});
                Files = dir(fullfile(FilePath,'*.csv'));
                if size(Files,1) > 1
                    error('Cannot have more than one .csv file in results time folder.  Please reduce the number of .csv files to 1 in each time folder containing results')
                end
                FullPath = fullfile(FilePath,Files.name);
                fid = fopen(FullPath, 'r');  %Open read only
                % Set headers and number of columns of data (24 columns currently)
                hdr{i} = textscan(fid, '%s',25,'delimiter',',');
                % Scan the numerical data (11 columns currently)
                bonds{i} = textscan(fid...
                    ,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'...
                    ,'delimiter',',','HeaderLines',1);
                fclose(fid);
                t1{i} = datevec(TimeFolds{i},'HHMM');
                idx(i) = etime(t1{i},t1{1})/60^2;
            end
            endtime = idx(i)+1;
            days = 1;
    end
end
end

