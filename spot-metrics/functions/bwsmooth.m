function [img] = bwsmooth(img, r)
%BWSMOOTH Smooth edges of bw areas as logical arrays

% based on https://www.mathworks.com/matlabcentral/answers/380687-how-to-smooth-rough-edges-along-a-binary-image

    %   img:    bw area as logical image
    %   r:      blur radius; will be rounded to odd for 2dconv

    if(nargin == 1)
        r = 7;
    else
        % round to nearest odd
        r = round((r-1)/2)*2+1;
    end
    
    kernel = ones(r) / r ^ 2;
    img = conv2(single(img), kernel, 'same');
    img = img > 0.5; % Rethreshold

end

