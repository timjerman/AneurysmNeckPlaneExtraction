function dist = distanceAcumulator(vertex, refPoints, refPointsRadius, scaleFactor)
% creates a volume of distances to closest vereces limited on the upper size
% by the estimated vessel endpoint radii

    volumeSize = round(scaleFactor*(max(vertex)-min(vertex))) + 1;

    refPoints = round(scaleFactor*(refPoints - repmat(min(vertex),size(refPoints,1),1))) + 1;
    refPointsRadius = refPointsRadius * scaleFactor - 2;
    refPointsRadius(refPointsRadius < 1.42) = 1.42;
    vertex = round(scaleFactor*(vertex - repmat(min(vertex),size(vertex,1),1))) + 1;
    
    dist = zeros(2*volumeSize);
    shift = round(volumeSize/2);
    
    refPoints = refPoints + repmat(shift,size(refPoints,1),1);    
    distRef1 = dist;
    distRef1(refPoints(1,1),refPoints(1,2),refPoints(1,3)) = 1;
    distRef1 = bwdist(distRef1);    
    distRef2 = dist;
    distRef2(refPoints(2,1),refPoints(2,2),refPoints(2,3)) = 1;
    distRef2 = bwdist(distRef2);
    
    rayEstRadiusLen = refPointsRadius(1) * distRef2 ./ (distRef1 + distRef2) + refPointsRadius(2) * distRef1 ./ (distRef1 + distRef2);
    
    idx = sub2ind(2*volumeSize,vertex(:,1)+shift(1),vertex(:,2)+shift(2),vertex(:,3)+shift(3));
    dist(idx) = 1;
    dist = bwdist(dist);
    dist(dist>rayEstRadiusLen) = rayEstRadiusLen(dist>rayEstRadiusLen);
    
    dist = dist(shift(1)+1:size(dist,1)-(volumeSize(1)-shift(1)),shift(2)+1:size(dist,2)-(volumeSize(2)-shift(2)),shift(3)+1:size(dist,3)-(volumeSize(3)-shift(3)));

    
    dist(~isfinite(dist)) = 0;
    
    