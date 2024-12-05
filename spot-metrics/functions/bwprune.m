function [roi] = bwprune(bw,c,r)
%BWPRUNE prunes thin attachments to area roughly centered at (c,r)

%   bw:     logical array containing area to be pruned
%   c,r:    column and row roughly corresponding to area center
%   out:    logical array with pruned area

    % initial selection in case the input image had multiple areas
    roi = bwselect(bw,c,r);
    
    % set up variables
    depth = 0;
    dist = bwdist(bwperim(roi));
    dist(~roi) = NaN;
    
    
    % until the roi is almost eroded
    while nnz(dist > depth) > 0.1*nnz(roi)
        % increase depth by one pixel
        depth = depth + 1;
        inner = dist > depth;
        
        % check if this has lead to an increase in components
        cc = bwconncomp(inner);
        if(cc.NumObjects > 1)
            
            % update the pointer to the area of interest if necessary
            if ~inner(c,r)
                % move it to the closest point that remains in the inner
                d = false(size(bw));
                d(c,r) = 1;
                d = bwdist(d);
                d(~inner) = NaN;
                [c,r] = find(d == min(min(d)));
            end
            
            % if so, select area of interest using pointer
            preselection = inner;
            difference = bwselect(preselection,c,r);
            % remove other components
            
            % dilate by depth
            difference = imdilate(...
                difference,strel('disk',depth,8));
            preselection = imdilate(...
                preselection,strel('disk',depth,8));
            
            % pruning step
            roi = bwselect(roi & ...
                ~bwconvhull( ...
                imdilate(preselection & ~difference, ...
                strel('disk',1,8))),c,r,8);
            
            % recompute distance transform
            depth = 0;
            dist = bwdist(bwperim(roi));
            dist(~roi) = NaN;
        end
    end

end
