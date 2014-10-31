function MechModuliPlot( EStarRaw, GStarRaw, nuStarRaw, EStarAr, GStarAr, nuStarAr,...
    YoungsRaw, LYoungsRaw, UYoungsRaw,...
    YoungsAr, LYoungsAr, UYoungsAr,...
    YoungsArTr, LYoungsArTr, UYoungsArTr,...
    ShearRaw, LShearRaw, UShearRaw,...
    ShearAr, LShearAr, UShearAr,...
    ShearArTr, LShearArTr, UShearArTr,...
    PoissonRaw, LPoissonRaw, UPoissonRaw,...
    PoissonAr, LPoissonAr, UPoissonAr,...
    PoissonArTr, LPoissonArTr, UPoissonArTr,...
    SRaw, SAr, SArTr, endtime )
% MechModuliPlot.m
% function MechModuliPlot( EStarRaw, GStarRaw, nuStarRaw, EStarAr, GStarAr, nuStarAr,...
%     YoungsRaw, LYoungsRaw, UYoungsRaw,...
%     YoungsAr, LYoungsAr, UYoungsAr,...
%     YoungsArTr, LYoungsArTr, UYoungsArTr,...
%     ShearRaw, LShearRaw, UShearRaw,...
%     ShearAr, LShearAr, UShearAr,...
%     ShearArTr, LShearArTr, UShearArTr,...
%     PoissonRaw, LPoissonRaw, UPoissonRaw,...
%     PoissonAr, LPoissonAr, UPoissonAr,...
%     PoissonArTr, LPoissonArTr, UPoissonArTr,...
%     SRaw, SAr, SArTr )
%
% This function plots the results of the function MechModuli.m.  It should
% be able to appropriately handle time series data.
%
% INPUTS:
%   EStar(---) - Scalar Young's modulus calculated from measurable
%   microstructural parameters.
%       *Raw - Results produced by the original contact tensor model.
%       *Ar - Results as produced by incorporating only the scalar shape
%       factor aspect ratio, Ar.
%
%   GStar(---) - Scalar shear modulus calculated from measurable
%   microstructural parameters.
%
%   nuStar(---) - Scalar Poisson's ratio calculated from measurable
%   microstructural parameters.
%
%   Youngs(---) - A structure of anisotropic Young's moduli along each of the 
%   principal directions of contact (i.e [E1,E2,E3]).
%       *Raw - Results as produced by the orginal contact tensor model.
%       *Ar - Results as produced by incorporating only the scalar shape
%       factor aspect ratio, Ar.
%       *ArTr - Results as produced by incorporating both the scalar shape
%       factor aspect ratio, Ar, and the tensor ratio, Tr.
%   
%   LYoungs(---)/UYoungs(---) - Lower and upper confidence intervals for
%   Youngs---.
%
%   Shear(---) - A structure of anisotropic shear moduli along aligning with 
%   the principal directions of contact (i.e. [G23,G13,G12]).
%
%   LShear(---)/UShear(---) - Lower and upper confidence intervals for
%   Shear(---).
%
%   Poisson(---) - A structure of anisotropic Poisson's ratio alignining with 
%   the principal directions of contact (i.e. [nu23,nu13,nu12]).
%
%   LPossion(---)/UPoisson(---) - Lower and upper confidence intervals for
%   Poisson(---).
%
%   S(---) - A 6x6 array that represents the anisotropic compliance matrix.
%   Its maximum deviation from isotropy is orthotropic symmetry (9 
%   independent coefficients) as restricted by the use of a 2nd order
%   contact tensor.
%
%   v1.0, October 30, 2014, 
% AUTHOR: David J. Walters; Montana State University

%% Raw
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Young''s Moduli, RAW',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsRaw.x./10^6,(YoungsRaw.x-LYoungsRaw.x)./10^6,...
    (UYoungsRaw.x-YoungsRaw.x)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsRaw.y./10^6,(YoungsRaw.y-LYoungsRaw.y)./10^6,...
    (UYoungsRaw.y-YoungsRaw.y)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsRaw.z./10^6,(YoungsRaw.z-LYoungsRaw.z)./10^6,...
    (UYoungsRaw.z-YoungsRaw.z)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarRaw./10^6,'k--')

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}','\it{E}\rm{_2}','\it{E}\rm{_3}','\it{E*}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli, RAW',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearRaw.xy./10^6,(ShearRaw.xy-LShearRaw.xy)./10^6,...
    (UShearRaw.xy-ShearRaw.xy)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearRaw.xz./10^6,(ShearRaw.xz-LShearRaw.xz)./10^6,...
    (UShearRaw.xz-ShearRaw.xz)./10^6,'bdiamond','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearRaw.yz./10^6,(ShearRaw.yz-LShearRaw.yz)./10^6,...
    (UShearRaw.yz-ShearRaw.yz)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar shear modulus G*
plot(idx,GStarRaw./10^6,'k--')

% Adjust format and appearance of shear modulus plot
legend('\it{G}\rm{_1_2}','\it{G}\rm{_1_3}','\it{G}\rm{_2_3}','\it{G*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poisson''s Ratio',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonRaw.xy,PoissonRaw.xy-LPoissonRaw.xy,...
    UPoissonRaw.xy-PoissonRaw.xy,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
legend('\it{\nu}\rm{_1_2}','\it{\nu*}');
grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonRaw.xz,PoissonRaw.xz-LPoissonRaw.xz,...
    UPoissonRaw.xz-PoissonRaw.xz,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
legend('\it{\nu}\rm{_1_3}','\it{\nu*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonRaw.yz,PoissonRaw.yz-LPoissonRaw.yz,...
    UPoissonRaw.yz-PoissonRaw.yz,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
legend('\it{\nu}\rm{_2_3}','\it{\nu*}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

%% Aspect Ratio, AR
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Young''s Moduli, Ar',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsAr.x./10^6,(YoungsAr.x-LYoungsAr.x)./10^6,...
    (UYoungsAr.x-YoungsAr.x)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsAr.y./10^6,(YoungsAr.y-LYoungsAr.y)./10^6,...
    (UYoungsAr.y-YoungsAr.y)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsAr.z./10^6,(YoungsAr.z-LYoungsAr.z)./10^6,...
    (UYoungsAr.z-YoungsAr.z)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'k--')

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}','\it{E}\rm{_2}','\it{E}\rm{_3}','\it{E*}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli, Ar',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearAr.xy./10^6,(ShearAr.xy-LShearAr.xy)./10^6,...
    (UShearAr.xy-ShearAr.xy)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearAr.xz./10^6,(ShearAr.xz-LShearAr.xz)./10^6,...
    (UShearAr.xz-ShearAr.xz)./10^6,'bdiamond','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearAr.yz./10^6,(ShearAr.yz-LShearAr.yz)./10^6,...
    (UShearAr.yz-ShearAr.yz)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar shear modulus G*
plot(idx,GStarAr./10^6,'k--')

% Adjust format and appearance of shear modulus plot
legend('\it{G}\rm{_1_2}','\it{G}\rm{_1_3}','\it{G}\rm{_2_3}','\it{G*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poisson''s Ratio',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonAr.xy,PoissonAr.xy-LPoissonAr.xy,...
    UPoissonAr.xy-PoissonAr.xy,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_2}','\it{\nu*}');
grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonAr.xz,PoissonAr.xz-LPoissonAr.xz,...
    UPoissonAr.xz-PoissonAr.xz,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_3}','\it{\nu*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonAr.yz,PoissonAr.yz-LPoissonAr.yz,...
    UPoissonAr.yz-PoissonAr.yz,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_2_3}','\it{\nu*}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

%% Aspect Ratio and Tensor Ratio, Ar, Tr
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Young''s Moduli, ArTr',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsArTr.x./10^6,(YoungsArTr.x-LYoungsArTr.x)./10^6,...
    (UYoungsArTr.x-YoungsArTr.x)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsArTr.y./10^6,(YoungsArTr.y-LYoungsArTr.y)./10^6,...
    (UYoungsArTr.y-YoungsArTr.y)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsArTr.z./10^6,(YoungsArTr.z-LYoungsArTr.z)./10^6,...
    (UYoungsArTr.z-YoungsArTr.z)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'k--')

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}','\it{E}\rm{_2}','\it{E}\rm{_3}','\it{E*}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli, ArTr',...
    'NumberTitle','off')
hold on

% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearArTr.xy./10^6,(ShearArTr.xy-LShearArTr.xy)./10^6,...
    (UShearArTr.xy-ShearArTr.xy)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearArTr.xz./10^6,(ShearArTr.xz-LShearArTr.xz)./10^6,...
    (UShearArTr.xz-ShearArTr.xz)./10^6,'bdiamond','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearArTr.yz./10^6,(ShearArTr.yz-LShearArTr.yz)./10^6,...
    (UShearArTr.yz-ShearArTr.yz)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar shear modulus G*
plot(idx,GStarAr./10^6,'k--')

% Adjust format and appearance of shear modulus plot
legend('\it{G}\rm{_1_2}','\it{G}\rm{_1_3}','\it{G}\rm{_2_3}','\it{G*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poisson''s Ratio',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonArTr.xy,PoissonArTr.xy-LPoissonArTr.xy,...
    UPoissonArTr.xy-PoissonArTr.xy,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_2}','\it{\nu*}');
grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonArTr.xz,PoissonArTr.xz-LPoissonArTr.xz,...
    UPoissonArTr.xz-PoissonArTr.xz,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_1_3}','\it{\nu*}');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonArTr.yz,PoissonArTr.yz-LPoissonArTr.yz,...
    UPoissonArTr.yz-PoissonArTr.yz,'ko','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarAr,'r--')
legend('\it{\nu}\rm{_2_3}','\it{\nu*}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

%% Plot All Together to Show Improvements of Each Model
% Plot tensor coefficients
% Initialize plot parameters
ebwidth = 0.1;  % Errorbar cap width
font = 'Palatino Linotype';
fsize = 11;     % Font Size
msize = 5;      % Marker Size

% Young's Modulus
% Initialize figure heading
figure('Name','3-D Young''s Moduli, All Together',...
    'NumberTitle','off')
hold on
% Plot Raw Values First
% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsRaw.x./10^6,(YoungsRaw.x-LYoungsRaw.x)./10^6,...
    (UYoungsRaw.x-YoungsRaw.x)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsRaw.y./10^6,(YoungsRaw.y-LYoungsRaw.y)./10^6,...
    (UYoungsRaw.y-YoungsRaw.y)./10^6,'rdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsRaw.z./10^6,(YoungsRaw.z-LYoungsRaw.z)./10^6,...
    (UYoungsRaw.z-YoungsRaw.z)./10^6,'ro','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarRaw./10^6,'r--')

% Next Plot Ar adjusted values
% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsAr.x./10^6,(YoungsAr.x-LYoungsAr.x)./10^6,...
    (UYoungsAr.x-YoungsAr.x)./10^6,'gsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsAr.y./10^6,(YoungsAr.y-LYoungsAr.y)./10^6,...
    (UYoungsAr.y-YoungsAr.y)./10^6,'gdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsAr.z./10^6,(YoungsAr.z-LYoungsAr.z)./10^6,...
    (UYoungsAr.z-YoungsAr.z)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'g--')

% Plot ArTr adjusted Values
% Plot value and confidence interval of Young's modulus E1
hE = errorbar(idx,YoungsArTr.x./10^6,(YoungsArTr.x-LYoungsArTr.x)./10^6,...
    (UYoungsArTr.x-YoungsArTr.x)./10^6,'bsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E2
hE = errorbar(idx,YoungsArTr.y./10^6,(YoungsArTr.y-LYoungsArTr.y)./10^6,...
    (UYoungsArTr.y-YoungsArTr.y)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of Young's modulus E3
hE = errorbar(idx,YoungsArTr.z./10^6,(YoungsArTr.z-LYoungsArTr.z)./10^6,...
    (UYoungsArTr.z-YoungsArTr.z)./10^6,'bo','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Adjust format and appearance of Young's modulus plot
legend('\it{E}\rm{_1}\rm{^R^a^w}','\it{E}\rm{_2}\rm{^R^a^w}','\it{E}\rm{_3}\rm{^R^a^w}',...
    '\it{E*}\rm{^R^a^w}',...
    '\it{E}\rm{_1}\rm{^A^r}','\it{E}\rm{_2}\rm{^A^r}','\it{E}\rm{_3}\rm{^A^r}',...
    '\it{E*}\rm{^A^r}',...
    '\it{E}\rm{_1}\rm{^A^r^T^r}','\it{E}\rm{_2}\rm{^A^r^T^r}','\it{E}\rm{_3}\rm{^A^r^T^r}',...
    'Location','NorthEast');
grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Young''s Modulus (MPa)');
% axis([-1 (endtime) 200 500])
xlim([-1 (endtime)])
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Shear Modulus
% Initialize figure heading
figure('Name','3-D Shear Moduli, All Together',...
    'NumberTitle','off')
hold on

% Plot Raw Values First
% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearRaw.xy./10^6,(ShearRaw.xy-LShearRaw.xy)./10^6,...
    (UShearRaw.xy-ShearRaw.xy)./10^6,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearRaw.xz./10^6,(ShearRaw.xz-LShearRaw.xz)./10^6,...
    (UShearRaw.xz-ShearRaw.xz)./10^6,'rdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearRaw.yz./10^6,(ShearRaw.yz-LShearRaw.yz)./10^6,...
    (UShearRaw.yz-ShearRaw.yz)./10^6,'ro','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarRaw./10^6,'r--')

% Next Plot Ar adjusted values
% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearAr.xy./10^6,(ShearAr.xy-LShearAr.xy)./10^6,...
    (UShearAr.xy-ShearAr.xy)./10^6,'gsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearAr.xz./10^6,(ShearAr.xz-LShearAr.xz)./10^6,...
    (UShearAr.xz-ShearAr.xz)./10^6,'gdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearAr.yz./10^6,(ShearAr.yz-LShearAr.yz)./10^6,...
    (UShearAr.yz-ShearAr.yz)./10^6,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value of scalar Young's Modulus E*
plot(idx,EStarAr./10^6,'g--')

% Plot ArTr adjusted Values
% Plot value and confidence interval of shear modulus G12
hE = errorbar(idx,ShearArTr.xy./10^6,(ShearArTr.xy-LShearArTr.xy)./10^6,...
    (UShearArTr.xy-ShearArTr.xy)./10^6,'bsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G13
hE = errorbar(idx,ShearArTr.xz./10^6,(ShearArTr.xz-LShearArTr.xz)./10^6,...
    (UShearArTr.xz-ShearArTr.xz)./10^6,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);

% Plot value and confidence interval of shear modulus G23
hE = errorbar(idx,ShearArTr.yz./10^6,(ShearArTr.yz-LShearArTr.yz)./10^6,...
    (UShearArTr.yz-ShearArTr.yz)./10^6,'bo','MarkerSize',msize);
errorbarwidth(hE,ebwidth);

% Adjust format and appearance of Young's modulus plot
legend('\it{G}\rm{_1}\rm{^R^a^w}','\it{G}\rm{_2}\rm{^R^a^w}','\it{G}\rm{_3}\rm{^R^a^w}',...
    '\it{G*}\rm{^R^a^w}',...
    '\it{G}\rm{_1}\rm{^A^r}','\it{G}\rm{_2}\rm{^A^r}','\it{G}\rm{_3}\rm{^A^r}',...
    '\it{G*}\rm{^A^r}',...
    '\it{G}\rm{_1}\rm{^A^r^T^r}','\it{G}\rm{_2}\rm{^A^r^T^r}','\it{G}\rm{_3}\rm{^A^r^T^r}',...
    'Location','NorthEast');

grid, box on
x1 = xlabel('Elapsed Time (hrs)');
y1 = ylabel('Shear Modulus (MPa)');
% axis([-1 (endtime) 65 95])
xlim([-1 (endtime)])
ylim('Auto')
set([y1 x1],'FontName',font,'FontSize',fsize)
set(gca,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)

% Poisson's Ratio
% Initialize figure heading
figure('Name','3-D Poisson''s Ratio, All Values',...
    'NumberTitle','off')

% Use subplots for better appearance
ax(1) = subplot(1,3,1);
hold on
% Plot value and confidence interval of Poisson's ratio nu12
hE = errorbar(idx,PoissonRaw.xy,PoissonRaw.xy-LPoissonRaw.xy,...
    UPoissonRaw.xy-PoissonRaw.xy,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonAr.xy,PoissonAr.xy-LPoissonAr.xy,...
    UPoissonAr.xy-PoissonAr.xy,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonArTr.xy,PoissonArTr.xy-LPoissonArTr.xy,...
    UPoissonArTr.xy-PoissonArTr.xy,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
plot(idx,nuStarAr,'g--')
title('\it{\nu}\rm{_1_2}')

grid, box on
y1 = ylabel('Poisson''s Ratio (-)');
axis([-1 (endtime) 0 0.5])


% Plot value and confidence interval of Poisson's ratio nu13
ax(2) = subplot(1,3,2);
hold on
hE = errorbar(idx,PoissonRaw.xz,PoissonRaw.xz-LPoissonRaw.xz,...
    UPoissonRaw.xz-PoissonRaw.xz,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonAr.xz,PoissonAr.xz-LPoissonAr.xz,...
    UPoissonAr.xz-PoissonAr.xz,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonArTr.xz,PoissonArTr.xz-LPoissonArTr.xz,...
    UPoissonArTr.xz-PoissonArTr.xz,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
plot(idx,nuStarAr,'g--')
title('\it{\nu}\rm{_1_3}')

grid, box on
x1 = xlabel('Elapsed Time (hrs)');
axis([-1 (endtime) 0 0.5])

% Plot value and confidence interval of Poisson's ratio nu23
ax(3) = subplot(1,3,3);
hold on
hE = errorbar(idx,PoissonRaw.yz,PoissonRaw.yz-LPoissonRaw.yz,...
    UPoissonRaw.yz-PoissonRaw.yz,'rsquare','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonAr.yz,PoissonAr.yz-LPoissonAr.yz,...
    UPoissonAr.yz-PoissonAr.yz,'go','MarkerSize',msize);
errorbarwidth(hE,ebwidth);
hE = errorbar(idx,PoissonArTr.yz,PoissonArTr.yz-LPoissonArTr.yz,...
    UPoissonArTr.yz-PoissonArTr.yz,'bdiamond','MarkerSize',msize+1);
errorbarwidth(hE,ebwidth);
plot(idx,nuStarRaw,'r--')
plot(idx,nuStarAr,'g--')
title('\it{\nu}\rm{_2_3}')

legend('\it{\nu}\rm{^R^a^w}',...
    '\it{\nu}\rm{^A^r}',...
    '\it{\nu}\rm{^A^r^T^r}',...
    '\it{\nu*}\rm{^R^a^w}','\it{\nu*}\rm{^A^r}');
grid, box on
axis([-1 (endtime) 0 0.5])

set([y1 x1],'FontName',font,'FontSize',fsize)
set(ax,'FontName',font,'FontSize',fsize,'XTick',0:3:endtime)
end

