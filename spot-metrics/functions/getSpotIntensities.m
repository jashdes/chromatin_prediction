function [intensities] = getSpotIntensities(img,pos,bg, nucleus)
%GETSPOTINTENSITIES Gets averaged, bg-subtracted, relative spot intensities

    %   img:    image from which intensities are to be extracted
    %   pos:    spot positions as n x 2 matrix
    %   bg:     logical array set for background pixels
    
    % initialize intensity array
    intensities = nan(size(pos,1),1);
    
    % convert image to double
    img = double(img);
    
    % perform background subtraction if a background has been passed
    if(nargin >= 3)
        
        % cut out foreground content
        background = img;
        background(~bg) = NaN;
        
        % calculate the median to estimate background noise
        noise = median(background, 'all', 'omitnan');
        
        % subtract background noise from remaining image
        img = img - noise;
        
    end
    
    % only consider nucleus of interest if map has ben passed
    if(nargin == 4)
        img(~nucleus) = NaN;
    end
    
    % relativize intensities
    
    img = img ./ mean(img,'all','omitnan');
    
    % approximate position to pixels
    pos = round(pos);
    
    sofar = false(size(img));
    
    for n = 1 : size(pos,1)
        
        % if the spot's coordinates are undefined
        if max(isnan(pos(n,:)))
            % its intensity is also undefined
            intensities(n) = NaN;
            continue;
        end
        
        % create a mask for this spot's position
        mask = false(size(img));
        mask(pos(n,1), pos(n,2)) = 1;
        
        % create area over which the intensities will be averaged
        mask = imdilate(mask,strel('square',3));
        
        sofar = sofar | mask;
        
        % calculate averages
        intensities(n) = mean(img(mask),'all');
    end
end

