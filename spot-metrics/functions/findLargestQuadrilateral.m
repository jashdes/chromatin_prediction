function[vertices] = findLargestQuadrilateral(pos)
%FINDLARGESTQUADRILATERAL Find vertices describing largest quadrilateral

    % pos:  [n x 2] grid tracking output where n is the number of spots in 
    %               the grid. Each column represents an axis
    % vertices:     rows of pos which were used as edge vertices

    % Get number of spots
    nspots = length(pos);
    % make sure it is sufficient
    assert(nspots >= 4);
    
    % Obtain possible spot combinations
    cmb = combnk(1:nspots,4);
    
    % Find vertices describing largest quadrilateral
    maxarea = 0;
    vertices = NaN(1,4);
    for c = 1 : length(cmb)
        spots = cmb(c,:);
        coordinates = [pos(spots(1),:); ...
                       pos(spots(2),:); ...
                       pos(spots(3),:); ...
                       pos(spots(4),:)];
        % If selection contains untracked spots, skip
        if max(isnan(coordinates(:))) == 1
            continue
        end
        x = coordinates(:,1);
        y = coordinates(:,2);
        chull = convhull(x,y);
        carea = polyarea(x(chull),y(chull));
        
        if carea > maxarea
            vertices = spots(chull);
            vertices = vertices(1:end-1);
            maxarea = carea;
            
        end
    end