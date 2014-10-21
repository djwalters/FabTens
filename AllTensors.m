function [PixSize,VolFrac,D,E,MeanStrucThick,StrucThickSD,F2,CMILRatio] = AllTensors(F2)
% TensorRoot.m
% This function is utilized for post processing computed fabric tensors.
% This reads the inputs from the calculated CONTACT FABRIC TENSOR and the
% MIL FABRIC TENSOR from CTAn.
%
% INPUTS:
%           F2: The second order contact tensor as calculated by the
%           function ContactTensor.m
%
% OUTPUTS:
%           PixSize: Outputs the pixel/voxel resolution of the CT Scan for
%           proper dimensional scaling.
%       
%           VolFrac: Ice volume fraction computed for the volume of
%           interest.
%
%           D: Eigen Values of the MIL Fabric Tensor reported from CTAn.
%
%           E: Eigen Vectors of the MIL Fabric Tensor reported from CTAn.
%
%           MeanStrucThick: The mean value of the structure thickness as
%           computed from CTAn.  This is produced by inscribing a sphere
%           along all voxels of the 3-D structures' skeleton.  The
%           structure thickness is the radius of that sphere.  In addition,
%           a histogram of the distribution of the computed structure
%           thicknesses is plotted.
%
%           StrucThickSD: Standard deviation of the structure thickness as
%           computed from CTAn.
%
%           F2: The 2nd order contact tensor as computed by ContactTensor.m
%
% VERSION: 1.0
% DATE: May 29, 2014
% AUTHOR: David J. Walters; Montana State University


%% Import data from CTAnData (external function)
[PixSize,VolFrac,D,E,MeanStrucThick,StrucThickHist,StrucThickSD]...
    = CTAnData;
%% Import data from Contact Tensor Segmentation Analysis
%[F2,~] = ContactTensor(1,PixSize); % NOTE: Second output of function is 4th order tensor if desired

%% Operate on Contact Fabric Tensor
% Generate an ellipsoid using the contact tensor to create a 3-D
% visualization of of the distribution of bonds in the analyzed sample.
% Generally assume contact tensor is in its principal orientation that is
% concurrently alligned with the global coordinate system.
% 
% To compare ellipsoids produced by the Contact Fabric Tensor, the MIL
% Fabric Tensor, and an isotropic sphere, all volumes are normalize to one.
VCON = (4/3)*pi*F2(1,1)*F2(2,2)*F2(3,3);
%Generates an ellipsoid with the radii scaled such that the volume equals 1
CRad1 = F2(1,1)/(VCON^(1/3));
CRad2 = F2(2,2)/(VCON^(1/3));
CRad3 = F2(3,3)/(VCON^(1/3));


[C1,C2,C3]=ellipsoid(0,0,0,CRad1,CRad2,CRad3,30);

%% Operate on the MIL Fabric Tensor
% Orient the MIL Fabric Tensor Ellipsoid with the Principal Axes
% D(1,1),D(2,2), and D(3,3) are the Eigenvalues (dimensions of the 
% ellipsoid axes)
% M1, M2, and M3 are the X, Y, and Z point clouds for producing a surface
% mesh.
% E is the EigenVector matrix as reported from CTAn.  Each column is the
% EigenVector associated with each EigenValue M1,M2,and M3
VMIL = (4/3)*pi*D(1)*D(2)*D(3);
%Generates an ellipsoid with the radii scaled such that the volume equals 1
MILRad1 = D(1)/(VMIL^(1/3));
MILRad2 = D(2)/(VMIL^(1/3));
MILRad3 = D(3)/(VMIL^(1/3));

[M1,M2,M3]=ellipsoid(0,0,0,MILRad1,MILRad2,MILRad3,30);
sz=size(M1);
for x=1:sz(1)
    for y=1:sz(2)
        A=[M1(x,y) M2(x,y) M3(x,y)]';
        A=E*A;
        M1(x,y)=A(1);M2(x,y)=A(2);M3(x,y)=A(3);
    end
end

%% Calculate Aspect ratio of between the MIL and Contact Ellipsoids
CMILRatio = zeros(3,3);
CMILRatio(1,1) = CRad1/MILRad1;
CMILRatio(2,2) = CRad2/MILRad2;
CMILRatio(3,3) = CRad3/MILRad3;

%% Compute volume of MIL Ellipsoid for equivalent isotropic sphere
rs = (3/4 / pi)^(1/3);

[S1,S2,S3] = sphere(30);
S1= S1*rs;
S2 = S2*rs;
S3 = S3*rs;

%% Plots
%Plot Ellipsoids of Contact Fabric Tensor, MIL Fabric Tensor, and Isotropic
%Sphere
figure('Name','3-D Tensor Comparisons','NumberTitle','off')
surf(C1,C2,C3,'FaceColor',[0 1 0],'EdgeColor',[0 0.5 0]);   %Green
hold on
surf(M1,M2,M3,'FaceColor',[0 0 1],'EdgeColor',[0 0 0.5]);   %Blue
surf(S1,S2,S3,'FaceColor',[1 0 0],'EdgeColor',[0.5 0 0]);   %Red
alpha(0.75)
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
legend('Contact Fabric Tensor','MIL Fabric Tensor','Isotropic Sphere'...
    ,'Location','best')

% Plot Distribution of Structure Thickness (Equivalent sphere radii) from
% CTAn to compare to Structure Thickness computed from sementation code.
figure('Name','Structure Thickness Distribution','NumberTitle','off')
bar(StrucThickHist(:,1),StrucThickHist(:,2),'hist');
xlabel('Pixels')
ylabel('Percent (%)')
YLimits = ylim;
hold on
line([MeanStrucThick MeanStrucThick], YLimits,'Color','r','LineWidth',2)
text(MeanStrucThick,YLimits(2)-(YLimits(2)*0.05)...
    ,'\leftarrow Mean Structure Thickness')

end