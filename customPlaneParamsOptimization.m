function bestParams = customPlaneParamsOptimization(vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir, lb0, ub0, initParams)

if nargin < 8
    initParams = [];
end

%%
if numel(initParams) == 0
    [px,py,pz] = meshgrid(linspace(lb0(1),ub0(1),3),linspace(lb0(2),ub0(2),9),linspace(lb0(3),ub0(3),50));
else
    [px,py,pz] = meshgrid(linspace(lb0(1),ub0(1),3),linspace(lb0(2),ub0(2),5),linspace(lb0(3),ub0(3),10));
end
paramSearch = [px(:) py(:) pz(:)];
paramSearch(paramSearch(:,1)==0,:) = [];
paramSearchUniqueZ = unique(pz(:));
paramSearch = [paramSearch;[zeros(numel(paramSearchUniqueZ),1) zeros(numel(paramSearchUniqueZ),1) paramSearchUniqueZ]];

if numel(initParams) > 0
    initParams(initParams < lb0) = lb0(initParams < lb0);
    initParams(initParams > ub0) = ub0(initParams > ub0);
    paramSearch = [paramSearch;initParams];
end

paramSearch = unique(paramSearch,'rows');

%%
solutions = zeros(size(paramSearch,1),1);

parfor i = 1 : size(paramSearch,1)
    cost = computeAneurysmNeckPlaneCostFuncMultipleCenterlines(paramSearch(i,:),vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir); 
    solutions(i,1) = cost;
end
%%
% [solutionsSorted, sortIdx] = sort(solutions,'ascend');
% tic
% [px,py,pz] = meshgrid(linspace(xBest(1)-10,xBest(1)+10,10),linspace(xBest(2)-25,xBest(2)+25,10),linspace(xBest(3)-0.05,xBest(3)+0.05,10));
% paramSearch = [px(:) py(:) pz(:)];
% solutions = zeros(size(paramSearch,1),1);
% for i = 1 : size(paramSearch,1)
%     cost = computeAneurysmNeckPlaneCostFuncMultipleCenterlines(paramSearch(i,:),vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir); 
%     solutions(i,1) = cost;
%     if costBest > cost
%         costBest = cost;
%         xBest = paramSearch(i,:);
%     end
% end
% toc

[solutionsSorted, sortIdx] = sort(solutions,'ascend');
paramSearchLast = paramSearch(sortIdx(1:5),:);

paramSearch = paramSearchLast;

for sI = 1 : 5
    lb = [max(lb0(1),paramSearchLast(sI,1)-10) max(lb0(2),paramSearchLast(sI,2)-25) max(lb0(3), paramSearchLast(sI,3) - 0.05)];
    ub = [min(ub0(1),paramSearchLast(sI,1)+10) min(ub0(2),paramSearchLast(sI,2)+25) min(max(lb0(3),paramSearchLast(sI,3)) + 0.05,ub0(3))];    
    
    [px,py,pz] = meshgrid(linspace(lb(1),ub(1),5),linspace(lb(2),ub(2),5),linspace(lb(3),ub(3),5));
    paramSearchTmp = [px(:) py(:) pz(:)];
    paramSearch = [paramSearch; paramSearchTmp];
end

%% test some solutions from genetic evolution 

paramSearchGenetic = [];
gaoptions = gaoptimset('Generations',1,'UseParallel','always','PopInitRange', [lb0; ub0],'Display','off');
for i = 1 : 10
    gasol = ga(@(x) computeAneurysmNeckPlaneCostFuncMultipleCenterlines(x,vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir),3,[],[],[],[],lb0,ub0,[],gaoptions);
    paramSearchGenetic = [paramSearchGenetic;gasol];
end

paramSearchGenetic = unique_tol(paramSearchGenetic, 1e-3,'rows');
paramSearch = [paramSearch;paramSearchGenetic];

%% remove unique parameters and remove multiple theta=0 angles
paramSearch = unique(paramSearch,'rows');
[~,uIdx,~] = unique(paramSearch(paramSearch(:,1)==0,[1,3]),'rows');
paramSearchZero = paramSearch(paramSearch(:,1)==0,:);
paramSearch = [paramSearchZero(uIdx,:); paramSearch(paramSearch(:,1)~=0,:)];

%%
solutions = zeros(size(paramSearch,1),1);
parfor i = 1 : size(paramSearch,1)
    cost = computeAneurysmNeckPlaneCostFuncMultipleCenterlines(paramSearch(i,:),vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir); 
    solutions(i,1) = cost;
end
%%

[solutionsSorted, sortIdx] = sort(solutions,'ascend');
paramSearchLast = paramSearch(sortIdx(1:10),:);
paramSearch = paramSearchLast;
for sI = 1 : 10
    lb = [max(lb0(1),paramSearchLast(sI,1)-5) max(lb0(2),paramSearchLast(sI,2)-10) max(lb0(3), paramSearchLast(sI,3) - 0.02)];
    ub = [min(ub0(1),paramSearchLast(sI,1)+5) min(ub0(2),paramSearchLast(sI,2)+10) min(max(lb0(3),paramSearchLast(sI,3)) + 0.02,ub0(3))];    
    
    [px,py,pz] = meshgrid(linspace(lb(1),ub(1),4),linspace(lb(2),ub(2),4),linspace(lb(3),ub(3),4));
    paramSearchTmp = [px(:) py(:) pz(:)];
    paramSearch = [paramSearch; paramSearchTmp];
end

% remove unique parameters and remove multiple theta=0 angles
paramSearch = unique(paramSearch,'rows');
[~,uIdx,~] = unique(paramSearch(paramSearch(:,1)==0,[1,3]),'rows');
paramSearchZero = paramSearch(paramSearch(:,1)==0,:);
paramSearch = [paramSearchZero(uIdx,:); paramSearch(paramSearch(:,1)~=0,:)];

%%
solutions = zeros(size(paramSearch,1),1);
parfor i = 1 : size(paramSearch,1)
    cost = computeAneurysmNeckPlaneCostFuncMultipleCenterlines(paramSearch(i,:),vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir); 
    solutions(i,1) = cost;
end
%%
[solutionsSorted, sortIdx] = sort(solutions,'ascend');
paramSearchLast = paramSearch(sortIdx(1:10),:);
paramSearch = paramSearchLast;
solutions = zeros(size(paramSearch,1),1);
lb = lb0;
ub = ub0;
parfor i = 1 : 10
    x0 = paramSearch(i,:);
%     lb = [0 -180 max(lb0(3), x0(3) - 0.05)];
%     ub = [45 180 min(max(lb0(3),x0(3)) + 0.05,1)];
%     options = psoptimset('UseParallel','always',...
    options = psoptimset(...
    'Display','off',...
                        'TolFun',1e-4,'TolX',1e-4,...
                        'InitialMeshSize', 100.0, ...
                        'CompleteSearch','on','SearchMethod',@searchlhs);
    [xg,fg] = patternsearch(@(x) computeAneurysmNeckPlaneCostFuncMultipleCenterlines(x,vertex,faces,vertexRefined,costMaps,costsNorm,costScaling,aneurysmPrincDir),x0,[],[],[],[],lb,ub,[],options);
    solutions(i,1) = fg;
    paramSearch(i,:) = xg;
end
%%
[bestCost, bestIdx] = min(solutions);
bestParams = paramSearch(bestIdx,:);

