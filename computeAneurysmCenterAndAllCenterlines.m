function [aneurysmCenter, centerlinesAll, aneurysmApex, meshInnerMask, radiusInterpFunc] = computeAneurysmCenterAndAllCenterlines(vertex, faces)
%% detect aneurysm center and all centerlines from a vascular mesh

accSize = 200;

%% resample mesh verteces

nRandSteps = 50;
ratioFaceSize = max(sqrt(sum((vertex(faces(:,1),:)-vertex(faces(:,2),:)).^2,2)))/mean(sqrt(sum((vertex(faces(:,1),:)-vertex(faces(:,2),:)).^2,2)));
nRandSteps = round(ratioFaceSize*nRandSteps);
nRandSteps(nRandSteps>200) = 200;
nRandSteps(nRandSteps<50) = 50;


rand_ratio = rand(nRandSteps,2);
while sum(rand_ratio(:,1)+rand_ratio(:,2)>=1)>0
    rand_ratio(rand_ratio(:,1)+rand_ratio(:,2)>=1,2) = rand(sum(rand_ratio(:,1)+rand_ratio(:,2)>=1),1);
end

vertexRefined = zeros(size(vertex,1) + (size(faces,1)*(nRandSteps)),3);
vertexRefined(1:size(vertex,1),:) = vertex;
for i = 1 : nRandSteps
    vertexRefined(size(vertex,1) + (i-1)*size(faces,1)+1:size(vertex,1) + (i)*size(faces,1),:) = vertex(faces(:,1),:) + rand_ratio(i,1)*(vertex(faces(:,2),:)-vertex(faces(:,1),:))+ rand_ratio(i,2)*(vertex(faces(:,3),:)-vertex(faces(:,1),:));
end

scaleFactor = round(accSize/max(abs(max(vertexRefined,[],1)-min(vertexRefined,[],1))));
volumeSize = round(scaleFactor*(max(vertexRefined)-min(vertexRefined))) + 1;

vertexRefined_ = round(scaleFactor*(vertexRefined - repmat(min(vertexRefined),size(vertexRefined,1),1))) + 1;
%% %%%%%%%%%%%%%%%%%%%%%   vessel endpoint detection
% find open sections of the vasculature
% asumption: vessels are open
%% detect edges - verteces where the vessel is cut

edges = [faces(:,[1,2]);faces(:,[1,3]);faces(:,[2,3])];
edges = sortrows(sort(edges,2));

[edgesuniq,ia,ic] = unique(edges,'rows');
[edgeshist] = hist(ic,1:max(ic));

edgesHole = edges(ia(edgeshist==1),:);
vertexHole = vertex(unique(edgesHole(:)),:);
%%
edgesHoleClustered = {[]};
currEIdx = edgesHole(1,:);
clusterIdx = 1;

while size(edgesHole,1) > 0    
    eRowIdx = any([repmat(edgesHole(:,1),1,numel(currEIdx)),repmat(edgesHole(:,2),1,numel(currEIdx))]==repmat(currEIdx,size(edgesHole,1),2),2);
    
    if sum(eRowIdx) == 0
        currEIdx = edgesHole(1,:);
        clusterIdx = clusterIdx + 1;
        edgesHoleClustered{clusterIdx} = [];
    else
        edgesHoleClustered{clusterIdx} = [edgesHoleClustered{clusterIdx}; edgesHole(eRowIdx,:)];
        currEIdx = edgesHole(eRowIdx,:);
        currEIdx = currEIdx(:)';
        edgesHole(eRowIdx,:) = [];
    end    
end
%%
vertexHoleClustered = {};
for i = 1 : numel(edgesHoleClustered)    
    vertexHoleClustered{i} = vertex(unique(edgesHoleClustered{i}(:)),:);
end
for i = 1 : numel(vertexHoleClustered)
    if numel(vertexHoleClustered{i}) < 5
        vertexHoleClustered(i) = [];
    end
end
display(['Number of vessel ends: ', num2str(numel(vertexHoleClustered))])

if numel(vertexHoleClustered) < 2
   error('Detected less than 2 vessel ends! Centerline cannot be constructed!'); 
end

%% compute and display normals for each vertex
[normal,~] = compute_normal(vertex,faces);
normal = -normal;

if size(normal,1) == 3
    normal = normal';
end

%% compute curvature
options.curvature_smoothing = 10;
options.verb = 1;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,~] = compute_curvature(vertex,faces,options);

% perform saturation on curvatures
CgaussSat = perform_saturation(Cgauss,1);
CgaussSat = rescale(max(min(CgaussSat,prctile(CgaussSat,95)),prctile(CgaussSat,5)),0,1);

CmeanSat = perform_saturation(Cmean,1);
CmeanSat = rescale(max(min(CmeanSat,prctile(CmeanSat,95)),prctile(CmeanSat,5)),0,1);

%figure,plot_fast_marching_mesh(vertex, faces, CgaussSat,{});    


%% compute the mesh inner mask -> volume inside the mesh

vertexDistMap = zeros(volumeSize);
idx = sub2ind(volumeSize,vertexRefined_(:,1),vertexRefined_(:,2),vertexRefined_(:,3));
vertexDistMap(idx) = 1;

%% a more uniform distribution of points
vertexHoleRefined = [];
numAdditionalVertecesPerHoleMin = 2000;
numAdditionalVertecesPerHoleMax = 10000;
numAdditionalVertecesPerHole = 2000;

for i = 1 : numel(vertexHoleClustered)
    
    holeSizeRatio = max(pdist(vertexHoleClustered{i}))/mean(sqrt(sum((vertex(faces(:,1),:)-vertex(faces(:,2),:)).^2,2)));    
    numAdditionalVertecesPerHole = round(holeSizeRatio * 4000 / 15);
    numAdditionalVertecesPerHole(numAdditionalVertecesPerHole<numAdditionalVertecesPerHoleMin) = numAdditionalVertecesPerHoleMin;
    numAdditionalVertecesPerHole(numAdditionalVertecesPerHole>numAdditionalVertecesPerHoleMax) = numAdditionalVertecesPerHoleMax;
    
    rand_ratio = rand(numAdditionalVertecesPerHole,2);
    while sum(rand_ratio(:,1)+rand_ratio(:,2)>=1)>0
        rand_ratio(rand_ratio(:,1)+rand_ratio(:,2)>=1,2) = rand(sum(rand_ratio(:,1)+rand_ratio(:,2)>=1),1);
    end
    
    vertexHoleRefinedTmp = zeros(numAdditionalVertecesPerHole, 3);
    for j = 1 : numAdditionalVertecesPerHole
%         vertexHoleRefinedTmp(j,:) = mean(vertexHoleClustered{i}(randperm(size(vertexHoleClustered{i},1),3)',:));
        randVert = randperm(size(vertexHoleClustered{i},1),3);
        vertexHoleRefinedTmp(j,:) = vertexHoleClustered{i}(randVert(1),:) + rand_ratio(j,1)*(vertexHoleClustered{i}(randVert(2),:)-vertexHoleClustered{i}(randVert(1),:))+ rand_ratio(j,2)*(vertexHoleClustered{i}(randVert(3),:)-vertexHoleClustered{i}(randVert(1),:));
    end
    vertexHoleRefined = [vertexHoleRefined;vertexHoleClustered{i}; vertexHoleRefinedTmp];
end
vertexHoleRefined_ = round(scaleFactor*(vertexHoleRefined - repmat(min(vertexRefined),size(vertexHoleRefined,1),1))) + 1;

%%
idx = sub2ind(volumeSize,vertexHoleRefined_(:,1),vertexHoleRefined_(:,2),vertexHoleRefined_(:,3));
vertexDistMap(idx) = 1;
vertexDistMap = bwdist(vertexDistMap);

%%
%%     detection of centerline endpoints
%%
%% weigtning function to prefer distant centerline points
nstart = 1;
clear options
options.W = ones(size(vertex,1),1);
options.nb_iter_max = Inf;
options.end_points = [];
pstarts = floor(rand(nstart,1)*size(vertex,1))+1;
[D,~,~] = perform_fast_marching_mesh(vertex', faces', pstarts, options);
%figure,plot_fast_marching_mesh(vertex,faces, D, [], options);
Df = D;
for i = 1 : 150
    pstarts = floor(rand(nstart,1)*size(vertex,1))+1;
    [D,~,~] = perform_fast_marching_mesh(vertex', faces', pstarts, options);
    Df = Df + D;
end
% figure,plot_fast_marching_mesh(vertex,faces, Df, [], options);

%% for each vertex send ray through the volume, detect the opposite vertex, compute the center point between both sides
% a possibility is to consider only those centers where the normals of both
% mesh sides have a small angle between them

rayCenters = [];
rayLengths = [];
rayCentersMinLengths = [];
rayCentersPropDist = [];
vertexOpposite = vertex * Inf;
rayCentersOrigIdx = [];

for vIdx = 1 : size(vertex,1)
%%
    orig = vertex(vIdx,:);
    dir = normal(vIdx,:);
    orig = orig + 0.01 * dir;
    
    if norm(dir) > 0

        [intersectRay,t,~,~,xcoor] = TriangleRayIntersection(orig, dir, vertex(faces(:,1),:), vertex(faces(:,2),:), vertex(faces(:,3),:),'lineType' , 'ray', 'border', 'inclusive');

        [ucoord,idx1,~] = unique(xcoor(intersectRay==1,:),'rows');
        tu = t(intersectRay==1,:);
        tu = tu(idx1);

    %% compute centers only if the normals have an opposite direction
        if numel(tu) > 0
            [~, minIdx] = min(tu);
            ucoord = ucoord(minIdx,:);            
            tu = tu(minIdx);
            vertexOpposite(vIdx,:) = ucoord;

            [~,minvIdx] = min(sqrt(sum((vertex-repmat(ucoord,size(vertex,1),1)).^2,2)));
            cosAngle = dot(normal(vIdx,:),-normal(minvIdx,:))/(norm(normal(vIdx,:))*norm(-normal(minvIdx,:)));        

            if cosAngle > cos(0.349066) %10 degrees = 0.174
                rayCenters = [rayCenters; orig + dir * tu/2];
                rayLengths = [rayLengths;tu];                
                [minVal,minIdx] = min(sqrt(sum((vertex-repmat(orig + dir * tu/2,size(vertex,1),1)).^2,2)));
                rayCentersMinLengths = [rayCentersMinLengths;minVal];
                rayCentersPropDist = [rayCentersPropDist; Df(minIdx)];
                
                rayCentersOrigIdx = [rayCentersOrigIdx; [vIdx minIdx]];
            end
        end
    end
end

%% retain centers where the detected radii is less then 1.5 times the 
% distance to the closest vertex

[~,searchDist] = knnsearch(vertex,rayCenters,'K',1);
rayCenters = rayCenters(0.5 * rayLengths < 1.5 * searchDist,:);
rayCentersPropDist = rayCentersPropDist(0.5 * rayLengths < 1.5 * searchDist);
rayCentersMinLengths = rayCentersMinLengths(0.5 * rayLengths < 1.5 * searchDist);
rayCentersOrigIdx = rayCentersOrigIdx(0.5 * rayLengths < 1.5 * searchDist,:);
rayLengths = rayLengths(0.5 * rayLengths < 1.5 * searchDist);

vertexIsParallel = zeros(size(vertex,1),1);
vertexIsParallel(unique(rayCentersOrigIdx(:))) = 1;
%% remove outliers and smooth

[IDX,searchDist] = knnsearch(rayCenters,rayCenters,'K',10);
rl = rayLengths(IDX);
rpd = rayCentersPropDist(IDX);
rml = rayCentersMinLengths(IDX);
distValid = searchDist(:,2:end) < repmat(0.5*rl(:,1),1,size(searchDist,2)-1);
mrl = sum(rl(:,2:end).*distValid,2) ./ sum(distValid,2);

gmrl = rl(:,1) > 1.5*mrl | searchDist(:,3) > 0.5*rl(:,1);
rayCenters = rayCenters(gmrl==0,:);
rl = rl(gmrl==0,:);
rpd = rpd(gmrl==0,:);
rml = rml(gmrl==0,:);
rayLengths = rayLengths(gmrl==0);
rayCentersPropDist = rayCentersPropDist(gmrl==0);
rayCentersMinLengths = rayCentersMinLengths(gmrl==0);
searchDist = searchDist(gmrl==0,:);

rpd(searchDist > repmat(rl(:,1)*0.5,1,size(rl,2))) = 0;
rml(searchDist > repmat(rl(:,1)*0.5,1,size(rl,2))) = 0;
rl(searchDist > repmat(rl(:,1)*0.5,1,size(rl,2))) = 0;

rayLengths = sum(rl,2)./sum(rl~=0,2);
rayCentersPropDist = sum(rpd,2)./sum(rpd~=0,2);
rayCentersMinLengths = sum(rml,2)./sum(rml~=0,2);
rayLengths(~isfinite(rayLengths)) = 0;
rayCentersPropDist(~isfinite(rayCentersPropDist)) = 0;
rayCentersMinLengths(~isfinite(rayCentersMinLengths)) = 0;
%%
rayCenters_ = round(scaleFactor*(rayCenters - repmat(min(vertexRefined),size(rayCenters,1),1))) + 1;
%%

%%
distFromHoles = min(pdist2(rayCenters,vertexHoleRefined),[],2);
distFromHoles = 1-distFromHoles/max(distFromHoles);
distFromHoles = distFromHoles * 1.15;
distFromHoles(distFromHoles>1) = 1;


%% based on the detected midpoint centers build a gaussian accumulator
lengthsScale = 0.5;
accAll = gaussianAcumulator(vertex, rayCenters, rayLengths*lengthsScale, scaleFactor,(distFromHoles).^2.*(2.*rayCentersMinLengths./rayLengths));

%% find maximums inside the accumulator and determine these points as real centers
accAll(~isfinite(accAll)) = 0;
accAll = imgaussian(accAll, 3);
%accAll = medfilt3(accAll,[5,5,5]);
mask = ones(3,3,3); mask(2,2,2) = 0;
accDil = imdilate(accAll, mask);
peaks = accAll > accDil;

[v1,v2,v3] = ind2sub(size(peaks),find(peaks>0));
endPoint_ = [v1,v2,v3];
endPoint = repmat(min(vertexRefined),size(endPoint_,1),1) + (endPoint_-1)/scaleFactor;

%%
% prepare interpolation function for radii extraction

radiusInterpFunc = scatteredInterpolant(rayCenters,0.5*rayLengths, 'nearest');

%%
%%

idx = sub2ind(volumeSize,endPoint_(:,1),endPoint_(:,2),endPoint_(:,3));
[~,meshInnerRefIdx] = max(vertexDistMap(idx));


%%
%%      mesh inner mask
%%
%%
meshInnerMask = ones(volumeSize);

meshInnerMask(vertexDistMap<1.42) = 0;
meshInnerMask(meshInnerMask>0) = 1;
meshInnerMask(~isfinite(meshInnerMask)) = 0;
meshInnerMask = single(~(meshInnerMask==0));

dt = DelaunayTri([vertexRefined_;vertexHoleRefined_]);  %# Create a Delaunay triangulation


[X,Y,Z] = ndgrid(1:size(meshInnerMask,1),1:size(meshInnerMask,2),1:size(meshInnerMask,3));   %# Create a mesh of coordinates for your volume
simplexIndex = pointLocation(dt,X(:),Y(:),Z(:));  %# Find index of simplex that
maskConvHull = ~isnan(simplexIndex);    %# Points outside the convex hull have a
maskConvHull = reshape(maskConvHull,size(X));   %# Reshape the mask to 101-by-101-by-101

meshInnerMask(maskConvHull==0) = 0;

CC = bwconncomp(meshInnerMask,6);
for i = 1 : numel(CC.PixelIdxList)
    meshInnerMask = meshInnerMask*0;
    meshInnerMask(CC.PixelIdxList{i}) = 1;
    if meshInnerMask(endPoint_(meshInnerRefIdx,1), endPoint_(meshInnerRefIdx,2), endPoint_(meshInnerRefIdx,3)) == 1
       break 
    end
end
meshInnerMask = imdilate(meshInnerMask,ones(3,3,3));
meshInnerMaskEroded = imerode(meshInnerMask,ones(5,5,5));
%%
%% based on the detected midpoint centers build a gaussian accumulator

accAll = accAll .* meshInnerMask;
mask = ones(3,3,3); mask(2,2,2) = 0;
accDil = imdilate(accAll, mask);
peaks = accAll > accDil;

[v1,v2,v3] = ind2sub(size(peaks),find(peaks>0));
endPoint_ = [v1,v2,v3];
endPoint = repmat(min(vertexRefined),size(endPoint_,1),1) + (endPoint_-1)/scaleFactor;
idx = sub2ind(volumeSize,endPoint_(:,1),endPoint_(:,2),endPoint_(:,3));

%%
%%      endpoint reduction
%%

endPoint(~meshInnerMask(idx),:) = [];
endPoint_(~meshInnerMask(idx),:) = [];

%% reduce endpoints to the number of vessel ends
endPointTmp = endPoint;
endPoint = zeros(numel(vertexHoleClustered),3);
for i = 1 : numel(vertexHoleClustered)
    edge2endpoindDist = pdist2(vertexHoleClustered{i},endPointTmp);
    [~,minIdx] = min(edge2endpoindDist,[],2);
    endPoint(i,:) = endPointTmp(mode(minIdx),:);
end
endPoint = unique(endPoint,'rows');
endPoint_ = round(scaleFactor*(endPoint - repmat(min(vertexRefined),size(endPoint,1),1))) + 1;
%%

if size(endPoint,1) < 2
   error('Detected less than 2 tentative centerline end points! Centerline cannot be constructed!'); 
end

%%
%% distance from endpoints -> vessel center is expected to be far away from the endpoints
%compute distance from the vessel end point - fastest path - to all points inside the mesh
% using fast marching (FM)

meshInnerMaskWeight = ones(volumeSize);

meshInnerMaskWeight(meshInnerMask==0) = 0.01;

W = double(meshInnerMaskWeight./max(meshInnerMaskWeight(:)));
W = rescale(W,1e-2,1);    

endPointFMDistanceMap = perform_fast_marching(W, endPoint_');
endPointFMDistanceMap = endPointFMDistanceMap .* meshInnerMask;
endPointFMDistanceMap = endPointFMDistanceMap/min(endPointFMDistanceMap(endPointFMDistanceMap(:)>0));

%%
%%     detection of all aneurysm centerlines
%%

%% compute distance from the aneurysm center point - fastest path - to all points inside the mesh
% % using fast marching (FM)

%% compute centerlines and remove those that are not the futhest in the selected direction

start_points_ = endPoint_;

currCount = 0;
centerlineData = [];
centerlinesAll = {};
pointsToRemove = [];

for startIdx = 1:size(start_points_,1)-1

    fm_start_point_ = start_points_(startIdx,:)';

    for endIdx = startIdx+1:size(start_points_,1)
     %%
        currCount = currCount + 1;
        
        fm_end_point_ = start_points_(endIdx,:)';
        
        %% compute distance map
        refPoints = endPoint([startIdx;endIdx],:);
        
        rayEstRadius = radiusInterpFunc(refPoints);
             
        accEstCent = distanceAcumulator(vertexRefined, refPoints, rayEstRadius, scaleFactor);

        
        accEstCent(meshInnerMask==0) = 0.01;
        
        accEstCentAll{currCount} = accEstCent;
        
        W = double(accEstCent./max(accEstCent(:)));
        W = rescale(W,1e-2,1);    
        weightedDistanceMap = perform_fast_marching(W, fm_start_point_);
        
%%
        options.method = 'continuous';
        minpath = compute_geodesic(weightedDistanceMap,fm_end_point_);
        if size(minpath,1) < size(minpath,2)
            minpath = minpath';
        end
        
        centerlineData = [centerlineData; [startIdx endIdx]];
        centerlinesAll{currCount} = (minpath-1) / scaleFactor + repmat(min(vertexRefined),size(minpath,1),1);       
        centerlinesAll{currCount} = unique(centerlinesAll{currCount},'rows','stable');        

    end
end

%%
pointsToRemove = unique(pointsToRemove);
rowToRem = unique([find(ismember(centerlineData(:,1),pointsToRemove));find(ismember(centerlineData(:,2),pointsToRemove))]);
centerlineData(rowToRem,:) = [];
centerlinesAll(rowToRem) = [];

%%
aneurysmCenterline = [];
for i = 1 : numel(centerlinesAll)
    aneurysmCenterline = [aneurysmCenterline; centerlinesAll{i}];
end
aneurysmCenterline = unique(aneurysmCenterline,'rows','stable');
aneurysmCenterline_ = round(scaleFactor*(aneurysmCenterline - repmat(min(vertexRefined),size(aneurysmCenterline,1),1))) + 1;
aneurysmCenterline_ = unique(aneurysmCenterline_,'rows','stable');

W = double(meshInnerMaskWeight./max(meshInnerMaskWeight(:)));
W = rescale(W,1e-2,1);

centerlineFMDistanceMap = perform_fast_marching(W, aneurysmCenterline_');
centerlineFMDistanceMap = centerlineFMDistanceMap .* meshInnerMask;
centerlineFMDistanceMap = centerlineFMDistanceMap/min(centerlineFMDistanceMap(centerlineFMDistanceMap(:)>0));

%%

%apexDistMap = endPointFMDistanceMap.*centerlineFMDistanceMap.*meshInnerMaskEroded;
apexDistMap = centerlineFMDistanceMap.*meshInnerMaskEroded;

[~, maxIdx] = max(apexDistMap(:));

[aneurysmApex(1,1),aneurysmApex(1,2),aneurysmApex(1,3)] = ind2sub(size(apexDistMap), maxIdx);
aneurysmApex = (aneurysmApex-1) / scaleFactor + repmat(min(vertexRefined),size(aneurysmApex,1),1);
aneurysmApex = vertex(knnsearch(vertex,aneurysmApex,'K',1),:);

%%
%%     detection of the aneurysm center
%%

%% find the closest point as the intersection between three randomly selected lines
% repeat many times, to randomly find such points
% -> aneurysm center

minIterations = 10000;
maxIterations = 100000;
minFoundPoints = 100;
minFoundPointsError = 5;
maxFoundPoints = 1000;

foundIntersections = [];
foundIntersectionsLengths = [];
foundIntersectionsTotalLength = [];
foundIntersectionsRadius = [];

for searchNeighborhood = [0.15, 0.5, 0.85]

    highCurvIdx = find(CmeanSat > 0.25 & ~vertexIsParallel);

    vertexSub = vertex(highCurvIdx,:);
    normalSub = normal(highCurvIdx,:);

    [IDX,~] = knnsearch(vertexSub,vertexSub,'K',round(searchNeighborhood*size(vertexSub,1)));
    randLims = round([0.05*size(IDX,2) size(IDX,2)]);

    iter = 0;
    while iter < maxIterations

        if (iter > minIterations && size(foundIntersections,1) > minFoundPoints) || size(foundIntersections,1) > maxFoundPoints
            break
        end

        iter = iter + 1;
        
        randIdx = randi(size(vertexSub,1));
        randIdx = [randIdx; IDX(randIdx, randi(randLims)); IDX(randIdx, randi(randLims))];
        vertexRand = [vertexSub(randIdx(1),:);vertexSub(randIdx(2),:);vertexSub(randIdx(3),:)];
        normalRand = [mean(normalSub(IDX(randIdx(1),1:6),:));mean(normalSub(IDX(randIdx(2),1:6),:));mean(normalSub(IDX(randIdx(3),1:6),:))];

        vertexRandLen = [sum((vertexRand(1,:) - vertexRand(2,:)).^2,2), sum((vertexRand(2,:) - vertexRand(3,:)).^2,2), sum((vertexRand(1,:) - vertexRand(3,:)).^2,2)];

        % all distances between points should be similar - prevents two very
        % close points to be considered
        if min(vertexRandLen)/max(vertexRandLen) > 0.6
            %%   
            normalRandAngles = dot(normalRand([1,2,3],:),normalRand([2,3,1],:),2);

            vertexRandDiff = vertexRand([2,3,1,3,1,2],:)-vertexRand([1,1,2,2,3,3],:);
            crossvec = cross(vertexRandDiff(1,:)/norm(vertexRandDiff(1,:)),vertexRandDiff(2,:)/norm(vertexRandDiff(2,:)));
            if all(cos(deg2rad(75)) < abs(dot(repmat(crossvec./norm(crossvec),3,1),normalRand,2)./arrayfun(@(x) norm(normalRand(x,:)),1:3)'))
            % angles between normals shouldn't be near 0 (perpendicular) or greater than 120
            % at least two normals should point in the same directions (two angles lower than 90)
                if all(normalRandAngles < cos(deg2rad(15))) && sum(normalRandAngles > cos(deg2rad(90))) > 1 && all(normalRandAngles > -0.5)
                    [P_intersect,distances,distancesToOrigin] = lineIntersect3D(vertexRand,normalRand);

                    normalBackProj = repmat(P_intersect,3,1)-vertexRand;
                    normalBackProj = normalBackProj ./ repmat(arrayfun(@(x) norm(normalBackProj(x,:)),1:3)',1,3);
                    angleBetweenOriginsDeg = rad2deg(acos(dot(vertexRandDiff([1,3,5],:),vertexRandDiff([2,4,6],:),2)./arrayfun(@(x) norm(vertexRandDiff(x,:)),[1,3,5])'./arrayfun(@(x) norm(vertexRandDiff(x,:)),[2,4,6])'));

                    % the angles between points shouldn't be too large and the
                    % intersections point shouldn't deviate too much from
                    % the initial lines
                    if all(distancesToOrigin>0) && all(angleBetweenOriginsDeg < 100) && all(cos(deg2rad(22.5)) < dot(normalRand,normalBackProj,2)./arrayfun(@(x) norm(normalRand(x,:)),1:3)')

                        distancesToOtherSide = zeros(3,1);
                        for nIdx = 1 : 3
                            orig = vertexRand(nIdx,:);
                            dir = normalRand(nIdx,:);
                            orig = orig + 0.01 * dir;
                            [intersectRay,t,~,~,xcoor] = TriangleRayIntersection(orig, dir, vertex(faces(:,1),:), vertex(faces(:,2),:), vertex(faces(:,3),:),'lineType' , 'ray', 'border', 'inclusive');

                            [~,idx1,~] = unique(xcoor(intersectRay==1,:),'rows');
                            tu = t(intersectRay==1,:);
                            tu = tu(idx1);
                            if numel(tu) > 0
                                %distancesToOtherSide(nIdx) = tu(1);
                                distancesToOtherSide(nIdx) = min(tu);
                            end
                        end    

                        % the intersection shouldn't be to close to meshe's other side
                        % and the distance between line to intersection
                        % compared to distances between points should be
                        % smaller
                        if all(distancesToOtherSide>0) && all(distancesToOrigin<0.6*distancesToOtherSide) && all(max(distances)/min(vertexRandLen) < 0.75)% && all((distances./distancesToOrigin)<0.25)           


                            if any(sum((P_intersect-repmat([1.775,-0.1691,1.499],size(P_intersect,1),1)).^2,2)<1e-2)
                               error('stop here') 
                            end                                

                                foundIntersections = [foundIntersections; P_intersect];
                                foundIntersectionsLengths = [foundIntersectionsLengths; mean(distancesToOrigin)];
                                foundIntersectionsTotalLength = [foundIntersectionsTotalLength; mean(distancesToOtherSide)];
                                foundIntersectionsRadius = [foundIntersectionsRadius; radiusInterpFunc(P_intersect)];        
                        end
                    end
                end
            end
        end
    end
    
    if size(foundIntersections,1) > minFoundPoints/4
       break 
    end
end
%%
display(['Center search: found intesections ', int2str(round(size(foundIntersections,1)) )])
if size(foundIntersections,1) < minFoundPoints
   warning(['Center search: number of found intesections ', int2str(round(size(foundIntersections,1)) ),' is bellow the recommended ',int2str(minFoundPoints),'. Solution: increase maximum number of iterations.']);
end
if size(foundIntersections,1) < minFoundPointsError
   error(['Center search: number of found intesections ', int2str(round(size(foundIntersections,1)) ),' is bellow the recommended ',int2str(minFoundPointsError),'. Solution: increase maximum number of iterations.']);
end

%% (new method) find the center of clusters by accumulating gaussians

%accCenter = gaussianAcumulator(vertex, foundIntersections, 2*foundIntersectionsLengths, scaleFactor);
% top: use estimated distance from intersection to mesh as radii
% bottom: use radii from interpolated radii map 
accCenter = gaussianAcumulator(vertex, foundIntersections, 2*foundIntersectionsRadius, scaleFactor);

mask = ones(3,3,3); mask(2,2,2) = 0;
accDil = imdilate(accCenter, mask);
peaks = accCenter > accDil & accCenter > max(accCenter(:))/10;

[v1,v2,v3] = ind2sub(size(peaks),find(peaks>0));
aneurysmCenterMaxPoint_ = [v1,v2,v3];
aneurysmCenterMaxPoint = repmat(min(vertexRefined),size(aneurysmCenterMaxPoint_,1),1) + (aneurysmCenterMaxPoint_-1)/scaleFactor;
% scatter3(aneurysmCenterMaxPoint(:,1),aneurysmCenterMaxPoint(:,2),aneurysmCenterMaxPoint(:,3),400,'filled')



%% take the tentative aneurysm center point that is furthest from the edpoints

idx = sub2ind(volumeSize,aneurysmCenterMaxPoint_(:,1),aneurysmCenterMaxPoint_(:,2),aneurysmCenterMaxPoint_(:,3));
[~, maxIdx] = max((apexDistMap(idx)./max(apexDistMap(:))));

aneurysmCenter = aneurysmCenterMaxPoint(maxIdx,:);
aneurysmCenter_ = round(scaleFactor*(aneurysmCenter - min(vertexRefined))) + 1;


%%
% figure; hold on;
% scatter3(vertex(:,1),vertex(:,2),vertex(:,3),20,'filled')
% scatter3(aneurysmCenterMaxPoint(:,1),aneurysmCenterMaxPoint(:,2),aneurysmCenterMaxPoint(:,3),200,maxPointVolRatio','filled')
% scatter3(aneurysmCenter(:,1),aneurysmCenter(:,2),aneurysmCenter(:,3),200,'filled')
