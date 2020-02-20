function [aneurysmCenterlinePoint, aneurysmPrincDir] = computeAneurysmPrincipalDirectionMultipleCenterlines(vertex, meshInnerMask, aneurysmCenter, aneurysmCenterlines, aneurysmApex)

computeToApex = 1;
if nargin < 5
    computeToApex = 0;
end

%%
aneurysmCenterline = [];
for i = 1 : numel(aneurysmCenterlines)
    aneurysmCenterline = [aneurysmCenterline; aneurysmCenterlines{i}];
end

%%
accSize = 200;
scaleFactor = round(accSize/max(abs(max(vertex,[],1)-min(vertex,[],1))));
volumeSize = size(meshInnerMask);

if computeToApex
    start_point = aneurysmApex;
else
    start_point = aneurysmCenter;
end
start_point = round(scaleFactor*(start_point - repmat(min(vertex),size(start_point,1),1))) + 1;
start_point = start_point';

accDist = distanceAcumulatorFromMesh(vertex, scaleFactor);

% accDist(accDist < 7) = 0;
accDist = meshInnerMask .* accDist;
W = double(accDist./max(accDist(:)));
W = rescale(W,1e-2,1);    
weightedDistanceMap = perform_fast_marching(W, start_point);

%
aneurysmCenterline_ = round(scaleFactor*(aneurysmCenterline - repmat(min(vertex),size(aneurysmCenterline,1),1))) + 1;
idx = sub2ind(volumeSize, aneurysmCenterline_(:,1),aneurysmCenterline_(:,2),aneurysmCenterline_(:,3));
distOfCenterlinePoints = weightedDistanceMap(idx);

[~,minIdx] = min(distOfCenterlinePoints);
aneurysmCenterlinePoint = aneurysmCenterline(minIdx,:);

%%

end_point = aneurysmCenterlinePoint;
end_point = round(scaleFactor*(end_point - repmat(min(vertex),size(end_point,1),1))) + 1;
end_point = end_point';

options.method = 'continuous';
minpath = compute_geodesic(weightedDistanceMap,end_point);
if size(minpath,1) < size(minpath,2)
    minpath = minpath'; 
end
minpathInACoord = minpath / scaleFactor + repmat(min(vertex),size(minpath,1),1);
minpathInACoord = unique(minpathInACoord,'rows','stable');


%% straight line in direction of the shortest path to the aneurysm center
% aneurysmPrincDir = minpathInACoord(2,:) - minpathInACoord(1,:);
% aneurysmPrincDir = aneurysmPrincDir ./ norm(aneurysmPrincDir);
%% straight line from aneurysmCenterlinePoint to aneurysm center
% aneurysmPrincDir = aneurysmCenter - aneurysmCenterlinePoint;
% aneurysmPrincDir = aneurysmPrincDir ./ norm(aneurysmPrincDir);

%% fastest path from principal direction to the centerline
aneurysmPrincDir = minpathInACoord;
numInterpPoints = 1000;
curveInt = linspace(0, 1, numInterpPoints);
curveInt = curveInt(1:end);
aneurysmPrincDir = interparc(curveInt,aneurysmPrincDir(:,1),aneurysmPrincDir(:,2),aneurysmPrincDir(:,3),'spline');

%% plot aneurysm principal direction
% 
% figure, hold on
% scatter3(vertex(:,1),vertex(:,2),vertex(:,3),20,'filled')
% scatter3(aneurysmCenterline(:,1),aneurysmCenterline(:,2),aneurysmCenterline(:,3),100,'filled')
% scatter3(aneurysmCenterlinePoint(1),aneurysmCenterlinePoint(2),aneurysmCenterlinePoint(3),300,'filled')
% scatter3(aneurysmCenter(1),aneurysmCenter(2),aneurysmCenter(3),300,'filled')
% scatter3(aneurysmApex(1),aneurysmApex(2),aneurysmApex(3),300,'filled')
% scatter3(aneurysmPrincDir(:,1),aneurysmPrincDir(:,2),aneurysmPrincDir(:,3),100,'filled')

