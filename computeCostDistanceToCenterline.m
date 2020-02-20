%% cost function 3: scaled distance from centerline
function costMap = computeCostDistanceToCenterline(vertex, aneurysmCenterlines, radiusInterpFunc)

aneurysmCenterline = [];
for i = 1 : numel(aneurysmCenterlines)
    aneurysmCenterline = [aneurysmCenterline; aneurysmCenterlines{i}];
end

aneurysmCenterline = unique(aneurysmCenterline,'rows','stable');

distToCent = zeros(size(vertex,1),1);

[~,vertexDist] = knnsearch(vertex,vertex,'K',2);
[shiftRangeX, shiftRangeY, shiftRangeZ] = meshgrid(-1:1, -1:1, -1:1);
shiftRange = normr([shiftRangeX(:) shiftRangeY(:) shiftRangeZ(:)]);
shiftRange(sum(shiftRange.^2,2)>1.05,:)=0;
shiftRange = [shiftRange * mean(vertexDist(:,2));shiftRange * mean(vertexDist(:,2)) * 0.5];

% slightly perturbe centerline points to avhieve smoothness of the mesh map
for i = 1 : size(shiftRange,1)  
%     randShift = (2*rand(size(aneurysmCenterline,1),3)-1)*shiftRange;
%     
%     [distToCentLidx,distToCentL] = knnsearch(aneurysmCenterline + randShift,vertex,'K',1);
%     [~,aneurysmCenterlineRadius] = knnsearch(vertex,aneurysmCenterline + randShift,'K',1);
%     distToCent = distToCent + (distToCentL-aneurysmCenterlineRadius(distToCentLidx));
    
    randShift = repmat(shiftRange(i,:),size(vertex,1),1);
    [distToCentLidx,distToCentL] = knnsearch(aneurysmCenterline,vertex + randShift,'K',1);
    
    distToCent = distToCent + (distToCentL-radiusInterpFunc(aneurysmCenterline(distToCentLidx,:)));
    
end
distToCent = distToCent / size(shiftRange,1)  ;
%%

% [distToCentLidx,~] = knnsearch(aneurysmCenterline,vertex,'K',1);
% [~,aneurysmCenterlineRadius] = knnsearch(vertex,aneurysmCenterline,'K',1);
% aneurysmCenterlineRadius = aneurysmCenterlineRadius(distToCentLidx);

[distToCentLidx,~] = knnsearch(aneurysmCenterline,vertex,'K',1);
aneurysmCenterlineRadius = radiusInterpFunc(aneurysmCenterline(distToCentLidx,:));

aneurysmCenterlineRadius = aneurysmCenterlineRadius ./ max(distToCent(:));
aneurysmCenterlineRadius = min(aneurysmCenterlineRadius,ones(size(aneurysmCenterlineRadius))*0.25);

distToCent = distToCent ./ max(distToCent(:));
distToCent(distToCent<=0) = 1e-5;
costMap = 1-exp(-(distToCent).^3./(2*aneurysmCenterlineRadius.^3));
% costMap = max(1-exp(-(distToCent).^3./(2*aneurysmCenterlineRadius.^3)), sqrt(distToCent));
%costMap = %1-exp(-(distToCent).^3/(2*0.25^3));

% figure,plot_fast_marching_mesh(vertex, faces,costMap,{}); 
