function [Tr,TrCi,MIL,MILCi] = TensorRatio(F2,F2Ci,MILVal,LMILVal,UMILVal,MILVec)
% TensorRatio.m
% Tr = TensorRatio(F2,MILVal,MILVec)
%
%   INPUTS:
%       F2(i,j,n): 2nd order contact tensor.  i and j reference the matrix
%       position within an individual tensor (i=1,2,3;j=1,2,3) and n
%       specifies the volume in chronological order with which the contact
%       tensor is associated.
%
%       MILVal: Eigen Values of the MIL Fabric Tensor reported from CTAn.
%
%       MILVec: Eigen Vectors of the MIL Fabric Tensor reported from CTAn.
%
%   OUTPUTS:
%       Tr: Tensor aspect ratio
%
%       MILRadVec: Vector containing MIL Values
%
% Version: 1.0 - October 28, 2014
%
% Version: 1.1 - Updated computation of the MIL tensor and reset the
% normaliztion parameters. October 29, 2015
%
% AUTHOR: David J. Walters; Montana State University

%% Ellipsoid Radii of Contact Tensor
% Generate an ellipsoid using the contact tensor to create a 3-D
% visualization of of the distribution of bonds in the analyzed sample.
% Generally assume contact tensor is in its principal orientation that is
% concurrently alligned with the global coordinate system.
%
% To compare ellipsoids produced by the Contact Fabric Tensor, the MIL
% Fabric Tensor, and an isotropic sphere, all volumes are normalize to one.
% Vcon = (4/3)*pi*F2(1,1)*F2(2,2)*F2(3,3);
%Generates an ellipsoid with the radii scaled such that the volume equals 1
% CRad1 = F2(1,1)/(Vcon^(1/3));
% CRad2 = F2(2,2)/(Vcon^(1/3));
% CRad3 = F2(3,3)/(Vcon^(1/3));
CRad1 = F2(1,1);
CRad2 = F2(2,2);
CRad3 = F2(3,3);

LCRad1 = F2Ci(1,1,1);
LCRad2 = F2Ci(2,2,1);
LCRad3 = F2Ci(3,3,1);

UCRad1 = F2Ci(1,1,2);
UCRad2 = F2Ci(2,2,2);
UCRad3 = F2Ci(3,3,2);

CRad = [CRad1; CRad2; CRad3];
LCRad = [LCRad1; LCRad2; LCRad3];
UCRad = [UCRad1; UCRad2; UCRad3];
%% Ellipsoid Radii of MIL Fabric Tensor

% Vmil = (4/3)*pi*MILVal(1)*MILVal(2)*MILVal(3);
%Generates an ellipsoid with the radii scaled such that the volume equals 1
% MILRad1 = MILVal(1)/(Vmil^(1/3));
% MILRad2 = MILVal(2)/(Vmil^(1/3));
% MILRad3 = MILVal(3)/(Vmil^(1/3));
MILRad1 = MILVal(1)/sum(MILVal);
MILRad2 = MILVal(2)/sum(MILVal);
MILRad3 = MILVal(3)/sum(MILVal);

LMILRad1 = LMILVal(1)/sum(LMILVal);
LMILRad2 = LMILVal(2)/sum(LMILVal);
LMILRad3 = LMILVal(3)/sum(LMILVal);

UMILRad1 = UMILVal(1)/sum(UMILVal);
UMILRad2 = UMILVal(2)/sum(UMILVal);
UMILRad3 = UMILVal(3)/sum(UMILVal);

% Adjust orientation of axes using the eigen vectors.  The Eigen Vector
% matrix should be 3x3 with the columns being individual vectors

MILRad1Vec = MILRad1 * MILVec(1,:);
MILRad2Vec = MILRad2 * MILVec(2,:);
MILRad3Vec = MILRad3 * MILVec(3,:);

LMILRad1Vec = LMILRad1 * MILVec(1,:);
LMILRad2Vec = LMILRad2 * MILVec(2,:);
LMILRad3Vec = LMILRad3 * MILVec(3,:);

UMILRad1Vec = UMILRad1 * MILVec(1,:);
UMILRad2Vec = UMILRad2 * MILVec(2,:);
UMILRad3Vec = UMILRad3 * MILVec(3,:);

MILRadVec = [MILRad1Vec;MILRad2Vec;MILRad3Vec];
LMILRadVec = [LMILRad1Vec;LMILRad2Vec;LMILRad3Vec];
UMILRadVec = [UMILRad1Vec;UMILRad2Vec;UMILRad3Vec];

MILSum = norm(MILRad1Vec) + norm(MILRad2Vec) + norm(MILRad3Vec);

% Tr(1) = CRad1/MILSum(1);
% Tr(2) = CRad2/MILSum(2);
% Tr(3) = CRad3/MILSum(3);

% Tr(1) = CRad1/norm(MILRad1Vec);
% Tr(2) = CRad2/norm(MILRad2Vec);
% Tr(3) = CRad3/norm(MILRad3Vec);

% Must match appropriate contact direction with oriented MIL direction
% since MIL Values 1, 2, and 3 don't necessarily point in the 1, 2, and 3
% directions.  The following loop figures out the orientation of each 
% principal MIL value and matches it with the correct principal contact
% value.
Index_prev = 0;
for i=1:3
    [~,Index] = max(abs(MILRadVec(i,:)));
    if i > 1
        if Index_prev(i) == Index || Index_prev(i-1) == Index
            if Index == 1
                [~,Index] = max(abs(MILRadVec(i,2:3)));
                Index = Index + 1;
            elseif Index == 2
                [~,Index] = max(abs(MILRadVec(i,1:2:3)));
                if Index == 2
                    Index = Index + 1;
                end
            elseif Index == 3
                [~,Index] = max(abs(MILRadVec(i,1:2)));
            end
        end
    end
%     Tr(Index) = norm(MILRadVec(i,:))/CRad(Index)
    Tr(Index) = CRad(Index)/norm(MILRadVec(i,:));
    LTr(Index) = LCRad(Index)/norm(LMILRadVec(i,:));
    UTr(Index) = UCRad(Index)/norm(UMILRadVec(i,:));
    
    MIL(Index,:) = MILRadVec(i,:);
    LMIL(Index,:) = LMILRadVec(i,:);
    UMIL(Index,:) = UMILRadVec(i,:);
    Index_prev(i+1) = Index;
end

TrCi(:,:,1) = LTr;
TrCi(:,:,2) = UTr;
MILCi(:,:,:,1) = LMIL;
MILCi(:,:,:,2) = UMIL;
% Tr(1) = norm(MILRad1Vec)/CRad1;
% Tr(2) = norm(MILRad2Vec)/CRad2;
% Tr(3) = norm(MILRad3Vec)/CRad3;