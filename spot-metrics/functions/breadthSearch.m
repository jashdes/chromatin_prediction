function [distancegraph] = breadthSearch(map,entrypoints)
%BREADTHSEARCH Perform breadth-first search on a logical matrix

% map: logical matrix representing the map to be searched
%       - accessible points are asserted
%       - inaccessible points are deasserted
% entrypoints: logical matrix of same shape as map
%       - points from which the search is started are asserted
% distancegraph: copy of map in which each pixel has a value equal to that
%                of its distance in pixels to the nearest entrypoint, or
%                NaN if no path exists to determine this distance

dim = size(map);

% only consider entry points that are accessible
entrypoints = entrypoints & map;

% collapse the entry points to indeces
entrypoints = find(entrypoints);

% initialize distance graph
distancegraph = NaN(dim);

% Use a Java PriorityQueue
import java.util.PriorityQueue;

% to import MapPoint successfully, make sure MATLAB can see MapPoint.class
% relevant documentation can be found at
% http://www.cs.yale.edu/homes/spielman/ECC/javaMatlab.html

% here, we assume it is located in the same folder as this function
javaaddpath(fileparts(mfilename('fullpath')));
import MapPoint.*;
q = PriorityQueue();

% Add all accessible entry points to the queue
for n = 1:numel(entrypoints)
    e = MapPoint(entrypoints(n), 0);
    add(q,e);
    distancegraph(entrypoints(n)) = 0;
end

% until the queue is empty
while ~(isEmpty(q))
    % pop an element off the queue
    current = poll(q);
    idx = current.index;
    dist = current.value;
    
    % define the list of its neighbors and their distance
    adjacency = [... 
        idx+1, 1; ... south
        idx-1, 1;... north
        idx + dim(1), 1; ... east
        idx - dim(1), 1; ... west
        idx - 1 - dim(1), sqrt(2); ... north-west
        idx - 1 + dim(1), sqrt(2); ... north-east
        idx + 1 - dim(1), sqrt(2); ... south-west
        idx + 1 + dim(1), sqrt(2); ... south-east
        ];
    
    % for each neighbor
    for n = 1:size(adjacency,1)
        % if it has not yet been visited and is accessible
        if isnan(distancegraph(adjacency(n,1))) && map(adjacency(n,1))
            % add it to the queue
            e = MapPoint(adjacency(n,1), dist + adjacency(n,2));
            add(q,e);
            % and mark it as visited
            distancegraph(adjacency(n,1)) = dist + adjacency(n,2);
        end
    end
    
end

end

