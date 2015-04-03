function MILPlot(idx,C1,C2,C3,M1,M2,M3,S1,S2,S3,...
    StrucThickHist,MeanStrucThick,StrucSepHist,MeanStrucSep,...
    PixSize,spatialLabel,VolFrac,endtime)
% MILPlot.m
% MILPlot(idx,C1,C2,C3,M1,M2,M3,S1,S2,S3,...
%     StrucThickHist,MeanStrucThick,PixSize,spatialLabel)
% This function is utilized for plotting data associated with the MIL
% fabric tensor and how it relates to the contact tensor
%
% INPUTS:
%
%   idx: Vector of elapsed time
%
%   C1,C2,C3: Matrices of X,Y,and Z data of contact tensor data produced
%   by the MATLAB Ellipsoid function for plotting with surf.
%
%   M1,M2,M3: Matrices of X,Y,and Z data of MIL tensor data produced
%   by the MATLAB Ellipsoid function for plotting with surf.
%
%   S1,S2,S3: Matrices of X,Y,and Z data of spheres produced by the MATLAB 
%   Ellipsoid function for plotting with surf.
%
%   StrucThickHist: Histogram data of structure thickness as computed from
%   CTAn.  This is produced by inscribing a sphere along all voxels of the 
%   3-D structures' skeleton.  The structure thickness is the radius of 
%   that sphere.  In addition, a histogram of the distribution of the 
%   computed structure thicknesses is plotted.
%
%   MeanStrucThick: The mean value of the structure thickness as
%   computed from CTAn.
%
%   StrucSepHist: Analogous to 'StrucThickHist' except operating on the
%   pore space of the snow microstructure
%
%   MeanStrucSep: The mean value of the structure seperation as computed
%   from CTAn.
%
%   PixSize: Outputs the pixel/voxel resolution of the CT Scan for
%   proper dimensional scaling.
%
%   spatialLabel: String variable set here to use in plotting functions
%   in other scripts and functions.
%
%   VolFrac: Volume fraction vector produced by CTAn
%
% Version: 1.0 - October 23, 2014
% AUTHOR: David J. Walters; Montana State University

% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size


%% Plot Ellipsoids of Contact Fabric Tensor, MIL Fabric Tensor, and Isotropic Sphere
for n = 1:length(idx)
str = sprintf('3-D Tensor Comparisons Time elapsed -- %2.0f hrs',idx(n));
figure('Name',str,'NumberTitle','off')
surf(C1{n},C2{n},C3{n},'FaceColor',[0 1 0],'EdgeColor',[0 0.5 0]);   %Green
hold on
surf(M1{n},M2{n},M3{n},'FaceColor',[0 0 1],'EdgeColor',[0 0 0.5]);   %Blue
surf(S1,S2,S3,'FaceColor',[1 0 0],'EdgeColor',[0.5 0 0]);   %Red
alpha(0.75)
axis equal
x1 = xlabel('X');
y1 = ylabel('Y');
z1 = zlabel('Z');
legend('Contact Fabric Tensor','MIL Fabric Tensor','Isotropic Sphere'...
    ,'Location','best')
set([z1 y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize)
end

%% Plot Distribution of Structure Thickness (Equivalent sphere radii) from
%CTAn to compare to Structure Thickness computed from sementation code.
for n = 1:length(idx)

str = sprintf('Structure Thickness Distribution Time elapsed -- %2.0f hrs',idx(n));
figure('Name',str,'NumberTitle','off')
bar(StrucThickHist{n}(:,1)*PixSize,StrucThickHist{n}(:,2),'hist');

x1 = xlabel(['Structure Thickness',spatialLabel]);
y1 = ylabel('Percent (%)');
YLimits = ylim;

hold on
line([MeanStrucThick(n)*PixSize MeanStrucThick(n)*PixSize],...
    YLimits,'Color','r','LineWidth',2)
text(MeanStrucThick(n)*PixSize,YLimits(2)-(YLimits(2)*0.05)...
    ,'\leftarrow      Mean Structure Thickness')
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize)
end

%% Plot Distribution of Structure Seperation (Equivalent sphere radii) from
%CTAn to compare to Structure Thickness computed from sementation code.
for n = 1:length(idx)

str = sprintf('Structure Seperation Distribution Time elapsed -- %2.0f hrs',idx(n));
figure('Name',str,'NumberTitle','off')
bar(StrucSepHist{n}(:,1)*PixSize,StrucSepHist{n}(:,2),'hist');

x1 = xlabel(['Structure Seperation',spatialLabel]);
y1 = ylabel('Percent (%)');
YLimits = ylim;

hold on
line([MeanStrucSep(n)*PixSize MeanStrucSep(n)*PixSize],...
    YLimits,'Color','r','LineWidth',2)
text(MeanStrucSep(n)*PixSize,YLimits(2)-(YLimits(2)*0.05)...
    ,'\leftarrow      Mean Structure Seperation')
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize)
end

%% Volume Fraction
figure('Name','Volume Fraction','NumberTitle','off')

% Plot volume fraction vs time
plot(idx,VolFrac,'ko','MarkerSize',msize);
grid
axis([-1 (endtime) 0 0.3])
y1 = ylabel('\phi: Ice Volume Fraction ');
x1 = xlabel('Elapsed Time (hrs)');

% Adjust final appearance of plots
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)
