function acc = gaussianAcumulator(vertex, rayCenters, rayLengths, scaleFactor, rayValue)
% creates a volume of accumulated  gaussians centered at rayCenters with
% sigma defined by rayLengths. The scale factor determines the size of the
% value, while rayValue applie per center gaussian scaling
%
%

    if nargin() < 5
        rayValue = ones(size(rayLengths));
    end

    volumeSize = round(scaleFactor*(max(vertex)-min(vertex))) + 1;

    rayCenters = round(scaleFactor*(rayCenters - repmat(min(vertex),size(rayCenters,1),1))) + 1;
    rayLengths = round(scaleFactor*(rayLengths));

    acc = zeros(2*volumeSize);
%     accN = acc;
    shift = round(volumeSize/2);

    for cIdx = 1 : size(rayCenters,1)    
        siz =  rayLengths(cIdx) * ones(1,3);
        sig = siz/(4*sqrt(2*log(2)));
        siz   = floor((siz-1)/2);%[1 1 1];%
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));    

        rayHBorders = [shift;shift] + [rayCenters(cIdx,:) - siz; rayCenters(cIdx,:) + siz];

        acc(rayHBorders(1,1):rayHBorders(2,1),rayHBorders(1,2):rayHBorders(2,2),rayHBorders(1,3):rayHBorders(2,3)) ...
            = acc(rayHBorders(1,1):rayHBorders(2,1),rayHBorders(1,2):rayHBorders(2,2),rayHBorders(1,3):rayHBorders(2,3)) ...
            + rayValue(cIdx) * h;
%         accN(rayHBorders(1,1):rayHBorders(2,1),rayHBorders(1,2):rayHBorders(2,2),rayHBorders(1,3):rayHBorders(2,3)) ...
%             = accN(rayHBorders(1,1):rayHBorders(2,1),rayHBorders(1,2):rayHBorders(2,2),rayHBorders(1,3):rayHBorders(2,3)) ...
%             + 1;
    end

%     accNmax = max(accN(:));
%     accN = accNmax*sqrt(accN./accNmax);
%     acc = acc ./ accN;    
    acc = acc(shift(1)+1:size(acc,1)-(volumeSize(1)-shift(1)),shift(2)+1:size(acc,2)-(volumeSize(2)-shift(2)),shift(3)+1:size(acc,3)-(volumeSize(3)-shift(3)));
    
    acc(~isfinite(acc)) = 0;
    %maxPoint = shiftCenters + (maxPoint-shift)/scaleFactor;
    
    
