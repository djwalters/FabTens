function [PixSize,VolFrac,D,E,MeanStrucThick,idx,F2,F2Ci,...
    F4,F4Ci,bondRad,grainRad,meanBondRad,meanGrainRad,coordNum,shapeFac,Tr]...
    = AllTensors(FileIn)
    
% AllTensors.m
% [PixSize,VolFrac,D,E,MeanStrucThick,idx,CMILRatio,F2,F2Ci,...
%     F4,F4Ci,bondRad,grainRad,meanBondRad,meanGrainRad,coordNum,shapeFac]...
%     = AllTensors(FileIn)
% This function is utilized for post processing computed fabric tensors.
% This reads the inputs from the calculated CONTACT FABRIC TENSOR and the
% MIL FABRIC TENSOR from CTAn.
%
% INPUTS:
%       FileIn: specifies whether a single file is used or multiple files
%       chronologically ordered. 1 for single file, 2 for multiple files
%
% OUTPUTS:
%       PixSize: Outputs the pixel/voxel resolution of the CT Scan for
%       proper dimensional scaling.
%       
%       VolFrac: Ice volume fraction computed for the volume of interest.
%
%       D: Eigen Values of the MIL Fabric Tensor reported from CTAn.
%
%       E: Eigen Vectors of the MIL Fabric Tensor reported from CTAn.
%
%       MeanStrucThick: The mean value of the structure thickness as
%       computed from CTAn.  This is produced by inscribing a sphere along
%       all voxels of the 3-D structures' skeleton.  The structure
%       thickness is the radius of that sphere.  In addition, a histogram
%       of the distribution of the computed structure thicknesses is
%       plotted.
%
%       F2(i,j,n): 2nd order contact tensor.  i and j reference the matrix
%       position within an individual tensor (i=1,2,3;j=1,2,3) and n
%       specifies the volume in chronological order with which the contact
%       tensor is associated.
%       
%       F2Ci(i,j,n,b): Upper and lower bonds of the 95% confidence interval
%       of 2nd order contact tensor.  i, j, and n are the same as above. 
%       b indicates the upper or lower bounds of the confidence interval.
%       1 indicates the lower bound, 2 indicates the upper bound.
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
%       Tr: Tensor ratio
%
% Version: 1.0 - May 29, 2014
% Version: 2.0 - October 23, 2014
%   - Added ability to view all time series data for both segmentation
%   contact tensor data and MIL tensor data produced by CTAn.
% Version: 2.1 - October 28, 2014
%   - Added function for calculating the Tensor Ratio.  Inclusion of
%   mechanical properties is the next task.
% AUTHOR: David J. Walters; Montana State University


%% Import data from CTAnData (external function)
[idx,PixSize,VolFrac,D,E,MeanStrucThick,StrucThickHist,RootPath]...
    = CTAnData(FileIn);
%% Import data from Contact Tensor Segmentation Analysis
[F2,F2Ci,F4,F4Ci,bondRad,grainRad,meanBondRad,meanGrainRad,...
    coordNum,shapeFac,idx,endtime,spatialLabel]...
    = ContactTensor(FileIn,1,PixSize,RootPath);
%% Operate on Contact Fabric Tensor
for n = 1:length(idx)
    idx(n)
    % Generate an ellipsoid using the contact tensor to create a 3-D
    % visualization of of the distribution of bonds in the analyzed sample.
    % Generally assume contact tensor is in its principal orientation that is
    % concurrently alligned with the global coordinate system.
    %
    % To compare ellipsoids produced by the Contact Fabric Tensor, the MIL
    % Fabric Tensor, and an isotropic sphere, all volumes are normalize to one.
    VCON(n) = (4/3)*pi*F2(1,1,n)*F2(2,2,n)*F2(3,3,n);
    %Generates an ellipsoid with the radii scaled such that the volume equals 1
    CRad1(n) = F2(1,1,n)/(VCON(n)^(1/3));
    CRad2(n) = F2(2,2,n)/(VCON(n)^(1/3));
    CRad3(n) = F2(3,3,n)/(VCON(n)^(1/3));
    
    
    [C1{n},C2{n},C3{n}]=ellipsoid(0,0,0,CRad1(n),CRad2(n),CRad3(n),30);
    
    %% Operate on the MIL Fabric Tensor
    % Orient the MIL Fabric Tensor Ellipsoid with the Principal Axes
    % D(1,1),D(2,2), and D(3,3) are the Eigenvalues (dimensions of the
    % ellipsoid axes)
    % M1, M2, and M3 are the X, Y, and Z point clouds for producing a surface
    % mesh.
    % E is the EigenVector matrix as reported from CTAn.  Each column is the
    % EigenVector associated with each EigenValue M1,M2,and M3
    VMIL(n) = (4/3)*pi*D(1,n)*D(2,n)*D(3,n);
    %Generates an ellipsoid with the radii scaled such that the volume equals 1
    MILRad1(n) = D(1,n)/(VMIL(n)^(1/3));
    MILRad2(n) = D(2,n)/(VMIL(n)^(1/3));
    MILRad3(n) = D(3,n)/(VMIL(n)^(1/3));
    
    [M1{n},M2{n},M3{n}]=ellipsoid(0,0,0,MILRad1(n),MILRad2(n),MILRad3(n),30);
    sz=size(M1{n});
    for x=1:sz(1)
        for y=1:sz(2)
            A=[M1{n}(x,y) M2{n}(x,y) M3{n}(x,y)]';
            A=E(:,:,n)*A;
            M1{n}(x,y)=A(1);M2{n}(x,y)=A(2);M3{n}(x,y)=A(3);
        end
    end
    
    

    
    %% Compute volume of MIL Ellipsoid for equivalent isotropic sphere
    rs = (3/4 / pi)^(1/3);
    
    [S1,S2,S3] = sphere(30);
    S1 = S1*rs;
    S2 = S2*rs;
    S3 = S3*rs;
    %% Calculate Aspect ratio of between the MIL and Contact Ellipsoids
    Tr(:,n) = TensorRatio(F2(:,:,n),D(:,n),E(:,:,n));
    
    %% Calculate Mechanical Properties
    % Be sure to include outputs for the original contact tensor model, the
    % outputs for the model modified just by the grain aspect ratio, and
    % the outputs for the model modified by both the grain aspect ratio and
    % tensor ratio.
    
    [EStarRaw(n), GStarRaw(n), nuStarRaw(n), EStarAr(n), GStarAr(n), nuStarAr(n),...
        YoungsRaw(n), LYoungsRaw(n), UYoungsRaw(n),...
        YoungsAr(n), LYoungsAr(n), UYoungsAr(n),...
        YoungsArTr(n), LYoungsArTr(n), UYoungsArTr(n),...
        ShearRaw(n), LShearRaw(n), UShearRaw(n),...
        ShearAr(n), LShearAr(n), UShearAr(n),...
        ShearArTr(n), LShearArTr(n), UShearArTr(n),...
        PoissonRaw(n), LPoissonRaw(n), UPoissonRaw(n),...
        PoissonAr(n), LPoissonAr(n), UPoissonAr(n),...
        PoissonArTr(n), LPoissonArTr(n), UPoissonArTr(n),...
        SRaw(:,:,n), SAr(:,:,n), SArTr(:,:,n)] = ...
        MechModuli(VolFrac(n), coordNum(n), meanBondRad(n), meanGrainRad(n), F2(:,:,n), F2Ci(:,:,n,:), shapeFac(n) , Tr(:,n));
    
end
%% Plots
%Change 0 back to idx
MILPlot(0,C1,C2,C3,M1,M2,M3,S1,S2,S3,...
    StrucThickHist,MeanStrucThick,PixSize,spatialLabel,VolFrac,endtime);
MechModuliPlot( EStarRaw, GStarRaw, nuStarRaw, EStarAr, GStarAr, nuStarAr,...
    YoungsRaw, LYoungsRaw, UYoungsRaw,...
    YoungsAr, LYoungsAr, UYoungsAr,...
    YoungsArTr, LYoungsArTr, UYoungsArTr,...
    ShearRaw, LShearRaw, UShearRaw,...
    ShearAr, LShearAr, UShearAr,...
    ShearArTr, LShearArTr, UShearArTr,...
    PoissonRaw, LPoissonRaw, UPoissonRaw,...
    PoissonAr, LPoissonAr, UPoissonAr,...
    PoissonArTr, LPoissonArTr, UPoissonArTr,...
    SRaw, SAr, SArTr, endtime, idx);

end