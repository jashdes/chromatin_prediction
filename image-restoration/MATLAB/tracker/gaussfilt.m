function [ imFilt ] = gaussfilt(im,sigma)
%Apply a Gaussian filter to a uniform time series
%   Inputs:im = image(numRows x numCols x numFrames) 
%       sigma = standard deviation of Gaussian filter to be applied.
%   Outputs: imFilt = image smoothed in time, i.e. over frames.
%   written by James Conder. Aug 22, 2013
%   convolution used for uniformly spaced t

% ADD SCALING
        [numRows, numCols, numFrames] = size(im);
        t=linspace(0, numFrames, numFrames+1);
        % n = length(z);  % number of data
a = 1/(sqrt(2*pi)*sigma);   % height of normalized Gaussian
sigma2 = sigma*sigma;

% check for uniform spacing
% if so, use convolution. If nonuniform in t, use numerical integration
uniform = false;
dt = diff(t);  % expecting dt = array of 1's for stack of images;
dt = dt(1);
ddiff = max(abs(diff(diff(t))));
if ddiff/dt < 1.e-4
    uniform = true;
end

if uniform
    % construct filter
    fsize = 3*sigma;    % for sigma = 3, fsize = 6 pixels
    xx = linspace(-fsize, fsize,2*fsize+1); % 2*fsize + 1 = 13,
        %  xx = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6] 
    filter = dt*a*exp(-0.5*((xx - mean(xx)).^2)/(sigma2));
    filter = filter/sum(filter);
    i = filter < dt*a*1.e-6;
    filter(i) = [];
    filterVector(:) = filter(1,:);
    sum(filterVector)
    
    % Scale the image intensities.
    for kk = 1:numFrames
        MM(kk) = max(max(im(:,:,kk)));
    end
    MMM = max(MM)
    for kk = 1:numFrames
        im_scaled(:,:,kk) = im(:,:,kk)./MMM;
    end
       
    for ii = 1:numRows
        for jj = 1:numCols
            imVector(:) = im_scaled(ii,jj,:);
            imVectorFiltScaled(:) = conv(imVector,filterVector,'same');  % MAIN OPERATION of gaussfilt
            imVectorFilt(:)       = imVectorFiltScaled.*MMM;             % unscale
            imFilt(ii,jj,:)= imVectorFilt(:);                     % vector -> matrix
        end
    end
%     for kk = 1:numFrames
%         imFilt(:,:,kk) = imFilt_scaled(ii,jj,kk).*MMM;
%     end
    
else
    fprintf('t is not uniformly spaced. Need to use gaussFilterOriginal')
       
end         % uniform vs non-uniform


end         % end of function

