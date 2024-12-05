function [] = imagepoint3Dplot(img,pos,data, radius)
%IMAGE3DPLOT Projects image onto X,Z plane with data on Y

    %   img:    image to project on X, Z plane
    %   pos:    n by 2 matrix with each row containing [X,Z] of datapoint
    %   data:   column vector containing the data to be projected
    %   radius: pixel radius to highlight data points
    
    
    % Set up the image plane
    [X,Z] = meshgrid(1:size(img,2),1:size(img,1));
    
    % add data to data plane
    Y= zeros(size(img));
    pos = round(pos);
    
    for k = 1:size(pos,1)
        if(max(isnan(pos(k,:))) == 0)
            Y(pos(k,1), pos(k,2)) = data(k);
        end
    end
    
    % set radius for highlighting
    if(nargin == 3)
        radius = 5;
    end
    
    if(radius ~= 0)
        % highlight data points
        Y = imdilate(Y,ones(radius,radius));
    end
    
    s = surf(X,Z,Y, 'CData',img(:,:,[1 1 1]),'FaceColor','texturemap');
    s.EdgeColor = 'none';
    set(gca,'Color','k');
    
    
end

