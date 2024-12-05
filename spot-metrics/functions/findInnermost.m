function [x,y, map] = findInnermost(map, scale)
%FINDINNERMOST Find shape center by minimizing average path distance to all border pixels

%   shape:  logical array containing the shape of object in question
%   scale:  factor by which to scale shape image; the closer to zero, the
%           faster but less accurate
%   x,y:    coordinates of innermost point
%   dstmap: AvgPathDstMap generated in computations; see getAvgPathDstMap.m

    if(nargin == 2)
        s = size(map);
        map = imresize(map, scale, 'nearest');
    end
    
    % get border of shape
    border = bwperim(map, 8);
    % generate average path distance map
    map = getAvgPathDstMap(map,border);
    % find minimum of this surface to get central point
    [x,y] = find(map == min(min(map)));
    x = x / scale;
    y = y / scale;
    if(nargin == 2)
        % interpolate map to its original size
        map = imresize(map, s);
    end
    
end

