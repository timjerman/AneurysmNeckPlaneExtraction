function costMap = computeCostIncurvation(vertex, aneurysmCenter, aneurysmCenterlines)
%%
aneurysmCenterline = [];
for i = 1 : numel(aneurysmCenterlines)
    aneurysmCenterline = [aneurysmCenterline; aneurysmCenterlines{i}];
end

aneurysmCenterline = unique(aneurysmCenterline,'rows','stable');

[~,vertexDist] = knnsearch(vertex,vertex,'K',2);
[shiftRangeX, shiftRangeY, shiftRangeZ] = meshgrid(-1:1, -1:1, -1:1);
shiftRange = normr([shiftRangeX(:) shiftRangeY(:) shiftRangeZ(:)]);
shiftRange(sum(shiftRange.^2,2)>1.05,:)=0;
shiftRange = [shiftRange * mean(vertexDist(:,2));shiftRange * mean(vertexDist(:,2)) * 0.5];

angleCentToCenter = zeros(size(vertex,1),1);
angleFactor = angleCentToCenter;

for i = 1 : size(shiftRange,1)   
    randShift = repmat(shiftRange(i,:),size(aneurysmCenterline,1),1);
    
    [distToCentLidx,~] = knnsearch(aneurysmCenterline + randShift,vertex,'K',1);
%     [distBetweenCentLidx,~] = knnsearch(aneurysmCenterline + randShift,aneurysmCenterline + randShift,'K',2);
    closestCentVertex = aneurysmCenterline(distToCentLidx,:);

    dirToACenter = repmat(aneurysmCenter,size(vertex,1),1) - vertex;
    dirToACenter = dirToACenter./repmat(sqrt(sum(dirToACenter.^2,2)),1,3);

    dirClosestCentVertex = vertex - closestCentVertex;
    dirClosestCentVertex = dirClosestCentVertex./repmat(sqrt(sum(dirClosestCentVertex.^2,2)),1,3);
    
    angleCentToCenter = angleCentToCenter + dot(dirToACenter, -dirClosestCentVertex,2);
    
%     crossDir = cross(dirToACenter,dirClosestCentVertex);
%     angleFactor = angleFactor + abs(dot(crossDir, normr(aneurysmCenterline(distToCentLidx,:)-aneurysmCenterline(distBetweenCentLidx(distToCentLidx,2),:)),2));
    
end

costMap = rescale((angleCentToCenter/size(shiftRange,1)+1)/2,0,1);
costMap = 1-exp(-costMap.^3/(2*0.5^3));
% costMap = costMap + asin(angleFactor/size(shiftRange,1))/pi*2;
% 
% numRand = 30;
% 
% aneurysmCenterline = [];
% for i = 1 : numel(aneurysmCenterlines)
%     aneurysmCenterline = [aneurysmCenterline; aneurysmCenterlines{i}];
% end
% 
% aneurysmCenterline = unique(aneurysmCenterline,'rows','stable');
% 
% [~,vertexDist] = knnsearch(vertex,vertex,'K',2);
% shiftRange = mean(vertexDist(:,2));
% 
% angleCentToCenter = zeros(size(vertex,1),1);
% 
% for i = 1 : numRand   
%     randShift = (2*rand(size(aneurysmCenterline,1),3)-1)*shiftRange;
%     
%     [distToCentLidx,~] = knnsearch(aneurysmCenterline + randShift,vertex,'K',1);
%     closestCentVertex = aneurysmCenterline(distToCentLidx,:);
% 
%     dirToACenter = repmat(aneurysmCenter,size(vertex,1),1) - vertex;
%     dirToACenter = dirToACenter./repmat(sqrt(sum(dirToACenter.^2,2)),1,3);
% 
%     dirClosestCentVertex = vertex - closestCentVertex;
%     dirClosestCentVertex = dirClosestCentVertex./repmat(sqrt(sum(dirClosestCentVertex.^2,2)),1,3);
%     
%     angleCentToCenter = angleCentToCenter + dot(dirToACenter, -dirClosestCentVertex,2);
% end
% 
% costMap = rescale((angleCentToCenter/numRand+1)/2,0,1);
% costMap = 1-exp(-costMap.^3/(2*0.5^3));