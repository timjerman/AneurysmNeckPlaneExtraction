%%
if exist('poolobj','var')
    delete(poolobj);
end
addpath('toolbox_fast_marching/')
addpath('toolbox_general/')
addpath('toolbox_graph/')

%% define case files - chien aneurysms
poolobj = parpool;

fileName = 'aneurysm.obj';

[vertex, faces] = readAneurysmMeshObj(fileName);
disp('Aneurysm mesh loaded!')

% merge same vertices
[vertex, ~, idxUnq] = unique(vertex, 'rows');
faces = idxUnq(faces);
faces = unique(faces,'rows');    

[aneurysmCenter, aneurysmCenterlines, aneurysmApex, meshInnerMask, radiusInterpFunc] = computeAneurysmCenterAndAllCenterlines(vertex, faces);
disp('Aneurysm centerline computed!')

[aneurysmCenterlinePoint, aneurysmPrincDir] = computeAneurysmPrincipalDirectionMultipleCenterlines(vertex, meshInnerMask, aneurysmCenter, aneurysmCenterlines, aneurysmApex);
disp('Aneurysm principal direction computed!')

costScaling = [1 1 1]; % (C_i, C_d, C_l) = (1 1 1)
[planeCurveCoord, planeParams, planeCurvePoints] = computeAneurysmPlaneParamsMultipleCenterlines(vertex,faces,aneurysmCenter,aneurysmPrincDir,aneurysmCenterlines,radiusInterpFunc,costScaling);

disp('Aneurysm plane optimized')    
disp(planeParams)

delete(poolobj);
%%
fig = figure;
plot_mesh(vertex, faces); 
alpha 0.8; 
material dull
shading('interp');
light;
lighting phong;
camlight('headlight');
hold on;       
scatter3(planeCurveCoord(:,1),planeCurveCoord(:,2),planeCurveCoord(:,3),150,'filled') 
