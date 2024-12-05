function [len, ang] = gridDistortion(pos)
%GRIDDISTORTION Measures the distortion two grid observations

    %   pos:    x grid observations with shape [x,n,2], where n is the
    %           number of spots in the grid. Each column represents an axis
    % length:   sum of absolute values of length changes of the largest
    %           quadrilateral with each vertex being a spot position
    %  angle:   sum of absolute values of angle changes of the the same
    %           quadrilateral

    % ensure number of spots is consistent across observations
    for n = 2:size(pos,1)
        assert(length(pos(n-1,:,:)) == length(pos(n,:,:)));
    end
    
    % remove grid corner point candidates that have interrupted tracking
    invalid = [];
    for k = 1 : size(pos,2)
        if max(isnan(pos(:,k,:)))
            invalid = [invalid, k];
        end
    end
    pos(:,invalid,:) = [];
            
    
    % Find spots that describe the largest quadrilateral
    vertices = findLargestQuadrilateral(squeeze(pos(1,:,:)));
    
    % Replace indices with position values
    p = pos(:,vertices,:);
    
    % Calculate vectors
    v = p(:,[2:end, 1],:) - p;
    % Close quadrilateral
    vc = v(:,[1:end,1], :);
    
    % Calculate lengths and angles
    lenv = NaN(size(v,1),4);
    angv = lenv;
    for n = 1 : size(v,1)
        for k = 1 : 4
            a = [squeeze(v(n,k,:));0];
            b = [squeeze(vc(n,k+1,:));0];
            lenv(n,k) = norm(a);
            angv(n,k) = atan2(norm(cross(a,b)), dot(a,b));
        end
    end
    
    % Calculate differences and sum the absolutes
    len = NaN(size(pos,1)-1,1);
    ang = len;
    for n = 1 : size(pos,1)-1
        len(n) = sum(abs(lenv(n+1,:) - lenv(n,:)),2);
        ang(n) = sum(abs(angv(n+1,:) - angv(n,:)),2);
    end
    
end
