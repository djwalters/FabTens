function TensorRatioPlot(idx,Tr,TrCi,endtime)
% This function reads the input from a seperate function, ContactTensor.m
% and plots the associated data for the production of useful figures and
% charts.  The inputs match the outputs of ContactTensor.m with the
% additional input of pixeldim to give real spatial resolution to the pixel
% based data.  The plots were originally inside of ContactTensor.m and have
% been moved to an external function for greater versatility.
%
% This function requires the use of errorbarwidth.m for appropriately
% sizing the error bars on the plots produced by this function.
%
% INPUTS
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
%       pixdim: Outputs the pixel/voxel resolution of the CT Scan for
%       proper dimensional scaling.
%
%       endtime: Variable listing the final elapsed time used in plotting
%       functions
%
%       spatialLabel: String variable set here to use in plotting functions
%       in other scripts and functions.
%
% Author: David J. Walters, Montana State University
% Version 1.0 - October 22, 2014
%   -Updated help file Oct. 23, 2014
%% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

%% Tensor Ratio
figure('Name','Tensor Ratio','NumberTitle','off')
% Plot value and confidence interval of contact tensor coefficient F11
hE = errorbar(idx,Tr(1,:),Tr(1,:)-TrCi(1,:,1),...
    TrCi(1,:,2)-Tr(1,:),'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of contact tensor coefficient F22
hold on
hE = errorbar(idx,Tr(2,:),Tr(2,:)-TrCi(2,:,1),...
    TrCi(2,:,2)-Tr(2,:),'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of contact tensor coefficient F33
hold on
hE = errorbar(idx,Tr(3,:),Tr(3,:)-TrCi(3,:,1),...
    TrCi(3,:,2)-Tr(3,:),'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Adjust format and appearance of contact tensor coefficient plot
grid
axis([-1 (endtime) 0.4 1.7])
y1 = ylabel('Tenosr Ratio Diagonals (-)');
x1 = xlabel('Elapsed Time (hrs)');
set([y1 x1],'FontName',font,'FontSize',fsize)
legend('\it{x}\rm{_1}','\it{x}\rm{_2}','\it{x}\rm{_3}',...
    'Location','best','Orientation','horizontal');
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% plot(idx,Tr(1,:),'rs','MarkerSize',msize);
% hold on
% plot(idx,Tr(2,:),'bd','MarkerSize',msize);
% plot(idx,Tr(3,:),'go','MarkerSize',msize);
% grid
% y1 = ylabel('Tensor Ratio (-)');
% x1 = xlabel('Elapsed Time (hrs)');
% set([y1 x1],'FontName',font,'FontSize',fsize)
% % axis([-1 (endtime) 1 3.5])
% xlim([-1 (endtime)])
% set(gca,'XTick',0:3:endtime)
% 
% legend('\it{x}\rm{_1}','\it{x}\rm{_2}','\it{x}\rm{_3}',...
%     'Location','Best','Orientation','horizontal');
% set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)