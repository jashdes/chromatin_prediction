function [p, normp, ratio] = getPositionMetrics(pos,center,nucleus)
%GETPOSITIONMETRICS Calculate position metrics for spot coordinates

%   pos:    matrix containing a two column [x,y] row for each spot
%   center: 1x2 matrix containing [x,y] of center
%   nucleus: logical array containing the shape of the nucleus
%   p:      row vector containing the distance to the nucleus border
%   normp:  row vector of p normalized by squareroot of area of nucleus
%   ratio:  row vector of c/(p+c) where c is distance to center

    % create empty distance map
    centerdistancemap = logical(zeros(size(nucleus)));
    
    % set relevant pixels
    centerdistancemap(center(1),center(2)) = 1; % only center pixel set
    borderdistancemap = bwperim(nucleus, 8);       % only border pixels set
    
    % calculate distance transform from set pixels and mask with nucleus
    centerdistancemap = bwdist(centerdistancemap) .* nucleus;
    borderdistancemap = bwdist(borderdistancemap) .* nucleus;
    
    % set pixels outside the nucleus to NaN as their distance is undefined
    centerdistancemap(centerdistancemap == 0) = NaN;
    borderdistancemap(borderdistancemap == 0) = NaN;
    
    % round positions
    pos = round(pos);
    
    % initialize vectors
    p = zeros(size(pos,1),1);
    c = zeros(size(pos,1),1);
    
    % for each spot
    for k = 1:size(pos,1)
        if(max(isnan(pos(k,:))) == 0)
            % get its distance from the nucleus border
            p(k) = borderdistancemap(pos(k,1),pos(k,2));
            % and the nucleus center
            c(k) = centerdistancemap(pos(k,1), pos(k,2));
        else
            p(k) = NaN;
            c(k) = NaN;
        end
    end
    
    % for all spots
    
    % calculate the normalized distance from the border
    normp = p./sqrt(bwarea(nucleus));
    
    % and the ratio
    ratio = c ./ (c + p);

end
