function [unitVecMat] = computeUnitVectorsOnCircle(nPts, typeFlag)

if (typeFlag == 1) % random draw
    thetaVec_rad = linspace(0,2*pi,nPts);
elseif (typeFlag == 2) % uniform
    thetaVec_rad = linspace(0,2*pi,nPts);
end

unitVecMat   = [cos(thetaVec_rad);sin(thetaVec_rad)];
