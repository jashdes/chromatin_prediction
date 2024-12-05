function [distancemap] = getAvgPathDstMap(map,entrypoints)
%FINDINNERMOST find point with least average distance from others

% map: logical matrix representing the map to be searched
%       - accessible points are asserted
%       - inaccessible points are deasserted
% entrypoints: logical matrix of same shape as map
%       - points from which the distance is calculated
% ip: point with lowest average distance from entrypoints

dim = size(map);

distancemap = zeros(dim);

% only consider entry points that are accessible
entrypoints = entrypoints & map;

% collapse the entry points to indeces
entrypoints = find(entrypoints);

% for each entry point
for n = 1 : numel(entrypoints)
    % create a logical image with only this entry point asserted
    e = zeros(dim);
    e(entrypoints(n)) = 1;
    e = logical(e);
    
    % and create a map of the distance of all points from this entrypoint
    summand = breadthSearch(map, logical(e));
    
    % replace NaN with 0
    summand(isnan(summand)) = 0;
    
    % and add it to the running sum
    distancemap = distancemap + summand;
end

distancemap(distancemap == 0) = NaN;

end

