function [costfunc, separateCosts, planeCurvePoints, meshPlaneIntersCoord, addedWeight] = computeAneurysmNeckPlaneCostFuncMultipleCenterlines(params,vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir,computeOnlyPlane)

if nargin < 9
    computeOnlyPlane = false;
    angleStep = 5;
else
    angleStep = 1;
end
%%
if numel(vertexRefined) == 0
   vertexRefined = vertex; 
end

%%
z = params(3);
phi = params(1);
theta = params(2);
%%
theta(theta>180) = theta(theta>180) - 360;
theta(theta<-180) = theta(theta<-180) + 360;

phi(phi>90) = phi(phi>90) - 180;
phi(phi<-90) = phi(phi<-90) + 180;

rotMat = azelaxes(theta,phi);

addedWeight = 0;

zIdx = round(z*(size(aneurysmPrincDir,1)-2)+1);
zIdx(zIdx<1) = 1;
zIdx(zIdx>(size(aneurysmPrincDir,1)-1)) = size(aneurysmPrincDir,1)-1;

aneurysmPrincDirOrig = aneurysmPrincDir(zIdx,:);
aneurysmPrincDirOrigDir = normr(aneurysmPrincDir(zIdx+1,:) - aneurysmPrincDir(zIdx,:));

dirCentToAneurysmPerp = normr(find_perp(aneurysmPrincDirOrigDir));
dirCentToAneurysmPerp2 = normr(cross(aneurysmPrincDirOrigDir,dirCentToAneurysmPerp));

B = [1,0,0;0,1,0;0,0,1];
A = [dirCentToAneurysmPerp;dirCentToAneurysmPerp2;aneurysmPrincDirOrigDir];
TM = B\A;

%%
angle = angleStep:angleStep: 360;
C = [cos(angle/180*pi);sin(angle/180*pi);0*angle];
C = (rotMat * C); 

dir = C'*TM;
dir = normr(dir);

meshPlaneIntersCoord = zeros(numel(angle),3);
meshPlaneIntersPoints = zeros(numel(angle),1);

for angleIdx =  1 : numel(angle)
    

    vertexRedux = 0.71 < dot(repmat(dir(angleIdx,:),size(vertex,1),1), normr(vertex-repmat(aneurysmPrincDirOrig,size(vertex,1),1)),2);
    vertexRedux = sort(find(vertexRedux == 1));
    % facesRedux = faces(any(reshape(builtin('_ismemberoneoutput',faces(:),vertexRedux),size(faces,1),3),2),:);
    facesRedux = faces(any(reshape(ismember(faces(:),vertexRedux),size(faces,1),3),2),:);
    
    
    if size(facesRedux,1) < 10
        facesRedux = faces;
    end
    
    [intersectRay,t,~,~,xcoor] = TriangleRayIntersection(aneurysmPrincDirOrig, dir(angleIdx,:), vertex(facesRedux(:,1),:), vertex(facesRedux(:,2),:), vertex(facesRedux(:,3),:),'lineType' , 'ray', 'border', 'inclusive','eps',1e-10);

%     [intersectRay,t,~,~,xcoor] = TriangleRayIntersection(aneurysmPrincDirOrig, dir(angleIdx,:), vertex(faces(:,1),:), vertex(faces(:,2),:), vertex(faces(:,3),:),'lineType' , 'ray', 'border', 'inclusive');

    [ucoord,idx1,~] = unique(xcoor(intersectRay==1,:),'rows');
    tu = t(intersectRay==1,:);
    tu = tu(idx1);

    if numel(tu) > 0
        %tu = tu(1);
        [~,minIdx] = min(tu);
        %tu = tu(minIdx);
        ucoord = ucoord(minIdx,:);
        meshPlaneIntersCoord(angleIdx,:) = ucoord;
        [~,minvIdx] = min(sqrt(sum((vertexRefined-repmat(ucoord,size(vertexRefined,1),1)).^2,2)));
        meshPlaneIntersPoints(angleIdx,1) = minvIdx;
    else
        addedWeight = addedWeight + 10;
    end    

end

idxRetain = all(meshPlaneIntersCoord~=0,2);
meshPlaneIntersPoints = meshPlaneIntersPoints(idxRetain);
meshPlaneIntersCoord = meshPlaneIntersCoord(idxRetain,:);

meshPlaneIntersCoord = unique(meshPlaneIntersCoord, 'rows','stable');

if size(meshPlaneIntersCoord,1) == 0 || size(meshPlaneIntersPoints,1) == 0
    costfunc = 1e10;
    separateCosts = [1e10 1e10 1e10];
    planeCurvePoints = [];
    meshPlaneIntersCoord = [];
    addedWeight = 1e10;
    return
end
%% cost functions

planeCurvePoints = vertexRefined(meshPlaneIntersPoints,:);

if computeOnlyPlane == false

    % length of the curve
%     len = (sqrt(sum(diff(vertex(meshPlaneIntersPoints,:)).^2,2)));
    len = (sqrt(sum(diff([meshPlaneIntersCoord;meshPlaneIntersCoord(1,:)]).^2,2)));
    costfuncLength = sum(len);    
    
    lenStretch = max(abs(diff([len(end);len])),abs(diff([len;len(1)])))./mean(len);
    
    if sum(lenStretch>2) > 0
        addedWeight = addedWeight + sum(lenStretch>2);
    end
    
    curveAngles = dot(normr(meshPlaneIntersCoord-[meshPlaneIntersCoord(end,:);meshPlaneIntersCoord(1:end-1,:)]),normr(meshPlaneIntersCoord-[meshPlaneIntersCoord(2:end,:);meshPlaneIntersCoord(1,:)]),2);
    if sum(curveAngles>=0) > 0
        addedWeight = addedWeight + sum(curveAngles>=0);
    end


    meshPlaneIntersPoints = unique(meshPlaneIntersPoints);
    
    costMapIncurvation = costMaps(:,1);
    costMapDistanceToCenterline = costMaps(:,2);
    
    costIncurvation = mean(costMapIncurvation(meshPlaneIntersPoints));
    costDistanceToCenterline = mean(costMapDistanceToCenterline(meshPlaneIntersPoints));

    separateCosts = [costIncurvation; costDistanceToCenterline; costfuncLength];

    Ndecimals = 4;
    f = 10.^Ndecimals;
    costfunc = costScaling(1) * round(f * costIncurvation) / f + costScaling(2) * round(f * costDistanceToCenterline)/ f +  costScaling(3) * round(f * costfuncLength/costsNorm(3))/f; 

    costfunc = costfunc + addedWeight;

    if ~isfinite(costfunc)
        costfunc = realmax;
    end
else
    costfunc = 0;
    separateCosts = 0;
end
