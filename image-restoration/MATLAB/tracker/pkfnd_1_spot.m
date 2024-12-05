function [CL nnn maxIntensity] = pkfnd_1_spot(im,nnn)
% replacement for pkfnd. Locates spot in small ROI for 7x7 V3.1 
% _2016_09_12  marks peak location in addition to profile line.
% INPUTS
    % im:  image 
    % nnn: figure number
% OUTPUT:  CL, 1 x 2 array giving location (x,y) for 1 spot.
%
% STEP 1.0 subtract baseline BL
% fprintf('pkfnd_1_spot');
% figure (nnn+1)
% subplot(1,3,1)
%     imshow(im,[]);
%     title('pkfnd-1-spot im');
    
BL = min(min(im));
im = im - BL;
% subplot(1,3,2)
%     imshow(im,[]);
%         title('pkfnd-1-spot im-BL')
% % find CM of spot
% fprintf('L22 pkfnd_1_spot');
maxIntensity= max(max(im));
mask  = (im>0.5*maxIntensity);
       
% apply mask to the image
    im_masked=im.*mask;
% subplot(1,3,3)
%    imshow(im_masked,[])
%       title('im-masked');
% pause (2)     
%           
CL = centerOfMass(im_masked); % Intensity-based center-of-mass location
     % CLl(2) gives x (cols); CL(1) gives y (rows). Subpixel.
end  % end of function
