function [planeCurveCoord, planeParams, planeCurvePoints] = computeAneurysmPlaneParamsMultipleCenterlines(vertex,faces,aneurysmCenter,aneurysmPrincDir,aneurysmCenterlines,radiusInterpFunc, costScaling, initParams)

if nargin < 8
    initParams = [];
end
if nargin < 7
    costScaling = [1 1 1];
end

%% refine mesh by subdividing each triangle into four smaller ones
% enables for higher resolution of cost function computation
% for ray to mesh intersection keep the original mesh because is way faster  
facesRefined = zeros(4*size(faces,1),3);
vertexRefined = ones(size(vertex,1)+3*size(faces,1),3)*nan;
vertexRefined(1:size(vertex,1),:) = vertex;
vertexRefinedCurrIdx = size(vertex,1);
for i = 1 : size(faces,1)
    newVertex = [mean(vertex(faces(i,[1,2]),:));mean(vertex(faces(i,[2,3]),:));mean(vertex(faces(i,[1,3]),:))];
    newVertexDist = pdist2(vertexRefined(1:vertexRefinedCurrIdx,:),newVertex);
    [minDist,newVertexIdx] = min(newVertexDist);
    vertexRefined(vertexRefinedCurrIdx+1:vertexRefinedCurrIdx+sum(minDist >= 1e-10),:) = newVertex(minDist >= 1e-10,:);
    newVertexIdx(minDist >= 1e-10) = vertexRefinedCurrIdx+1 : vertexRefinedCurrIdx+sum(minDist >= 1e-10);
    vertexRefinedCurrIdx = vertexRefinedCurrIdx + sum(minDist >= 1e-10);

    facesRefined((i-1)*4+1,:) = [faces(i,1) newVertexIdx(1) newVertexIdx(3)];
    facesRefined((i-1)*4+2,:) = [newVertexIdx(1) faces(i,2) newVertexIdx(2)];
    facesRefined((i-1)*4+3,:) = [newVertexIdx(3) newVertexIdx(1) newVertexIdx(2)];             
    facesRefined((i-1)*4+4,:) = [newVertexIdx(3) newVertexIdx(2) faces(i,3)];  
end
vertexRefined(isnan(vertexRefined(:,1)),:) = [];

%% cost function 1: incurvation - angle between vertex-centerline and centerline-aneurysm center vectors
costMapIncurvationRefined = computeCostIncurvation(vertexRefined, aneurysmCenter, aneurysmCenterlines);
% figure,plot_fast_marching_mesh(vertex, faces,costMapIncurvation,{});
%% cost function 2: scaled distance from centerline
costMapDistanceToCenterlineRefined = computeCostDistanceToCenterline(vertexRefined, aneurysmCenterlines, radiusInterpFunc);
% figure,plot_fast_marching_mesh(vertex, faces,costMapDistanceToCenterline,{});
%%
costMaps = cat(2,costMapIncurvationRefined,costMapDistanceToCenterlineRefined);
%%
z0 = 0.5;
[~, z0] = min(pdist2(aneurysmCenter, aneurysmPrincDir));
z0 = z0 / size(aneurysmPrincDir,1);
%%

[~, costsNorm,~,~,addedWeight] = computeAneurysmNeckPlaneCostFuncMultipleCenterlines([0 0 z0],vertex,faces,vertexRefined,costMaps,[1 1 1],costScaling,aneurysmPrincDir);
if addedWeight > 0
   for zStep = 0.1 : 0.1 : 0.9
       [ ~, costsNorm,~,~,addedWeight] = computeAneurysmNeckPlaneCostFuncMultipleCenterlines([0 0 zStep],vertex,faces,vertexRefined,costMaps,[1 1 1],costScaling,aneurysmPrincDir);
       if addedWeight == 0
           break
       end
   end
end
if addedWeight > 0
    costsNorm = [1,1,2*pi*radiusInterpFunc(aneurysmCenter)];
end
costsNorm = [1,1,costsNorm(3)];

%%

centerlineBoundLower = 0.25 * radiusInterpFunc(aneurysmPrincDir(1,:))/sum(sqrt(sum(diff(aneurysmPrincDir).^2,2)));
centerlineBoundUpper = 0.7;
if centerlineBoundLower > centerlineBoundUpper
    centerlineBoundLower = 0.6;
    centerlineBoundUpper = 0.8;
end
lb0 = [0,-180,centerlineBoundLower];
ub0 = [45,180,centerlineBoundUpper];
tic
planeParams = customPlaneParamsOptimization(vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir, lb0, ub0,initParams);
toc
[cost,costVals,planeCurvePoints,planeCurveCoord, addedWeight] = computeAneurysmNeckPlaneCostFuncMultipleCenterlines(planeParams,vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir);
