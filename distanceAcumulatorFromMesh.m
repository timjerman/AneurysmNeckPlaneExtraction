function dist = distanceAcumulatorFromMesh(vertex, scaleFactor)
% creates a volume of distances to closest vereces limited on the upper size
% by the estimated vessel endpoint radii

    volumeSize = round(scaleFactor*(max(vertex)-min(vertex))) + 1;
    vertex = round(scaleFactor*(vertex - repmat(min(vertex),size(vertex,1),1))) + 1;
    
    acc = zeros(2*volumeSize);
    shift = round(volumeSize/2);
        
    idx = sub2ind(2*volumeSize,vertex(:,1)+shift(1),vertex(:,2)+shift(2),vertex(:,3)+shift(3));
    acc(idx) = 1;
    dist = bwdist(acc);
    
    dist = dist(shift(1)+1:size(dist,1)-(volumeSize(1)-shift(1)),shift(2)+1:size(dist,2)-(volumeSize(2)-shift(2)),shift(3)+1:size(dist,3)-(volumeSize(3)-shift(3)));
 
    
    dist(~isfinite(dist)) = 0;
    
    