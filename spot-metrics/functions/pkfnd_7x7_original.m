function [out goodSpot] = pkfnd_7x7_original(im,th,sz,rr,goodSpot)
% finds local maxima to +/- 1 pixel precision.   
%  written by Eric Dufresne
%  Reworked by gholzwarth to assign rr systematically, using masked template  11/18/2016 by gholz

% INPUTS:
% im: image
% th: threshold minimum brightness of a pixel that might be a local maxima. 
% sz:  if your data is noisy, set this optional keyword to a value slightly
    % larger than the diameter of your blob. If multiple peaks are found
    % within a radius of sz/2 then the code will keep only the brightest.
    % Also gets rid of all peaks within sz of boundary
% rr = spotNum
% goodSpot: 49 x 5 matrix with 1's for good spots, 0's for bad spots
%           col1 = spotID; col2

% Imax_im_filtered is max(max(im_filtered(:,:,1)
% OUTPUT:  array containing x,y coordinates of local maximum rr 
%           out(rr,1) is the x-coordinate (column) of spot rr
%           out(rr,2) is the y-coordinate (rows) of spot rr
%           goodSpot modified if spot is too weak.
% CREATED: Eric R. Dufresne, Yale University, Feb 4 2005

% Modified by gholzwarth. Expects input images which are masked ROIs
    % with only 1 peak. 

%find all the pixels above threshold
%im=im./max(max(im));

ind=find(im > th);  % find indices of pixels at which im>threshold
[nr,nc]=size(im);

n=length(ind);
if n==0
    out=[];
    display('nothing above threshold');
    return;
end

mx=[];  % mx declared as array. Why did programmer call it "mx" when it has both rows and cols?

%convert index from find(L26) to row and column
rc=[mod(ind,nr),floor(ind/nr)+1]; % mod = modulus
for i=1:n
    r=rc(i,1);c=rc(i,2);
    %check each pixel above threshold to see if it is brighter than its neighbors
    if r>1 & r<nr & c>1 & c<nc
        if im(r,c)>=im(r-1,c-1) & im(r,c)>=im(r,c-1) & im(r,c)>=im(r+1,c-1) & ...
         im(r,c)>=im(r-1,c)  & im(r,c)>=im(r+1,c) &   ...
         im(r,c)>=im(r-1,c+1) & im(r,c)>=im(r,c+1) & im(r,c)>=im(r+1,c+1)
        mx=[mx,[r,c]']; % NOTE transpose. What does this line do?
        end
    end
end
mx=mx';  % NOTE second transpose.  Does code writer know what he's doing?

[npks,crap]=size(mx);  %npks is the number of peaks.

%if size is specified, then get rid of pks within size of boundary
if nargin==3 & npks>0
   %throw out all pks within sz of boundary;
    ind=find(mx(:,1)>sz & mx(:,1)<(nr-sz) & mx(:,2)>sz & mx(:,2)<(nc-sz));
    mx=mx(ind,:); 
end

%prevent from finding peaks within size of each other
[npks,crap]=size(mx);
if npks > 1 
    %CREATE AN IMAGE WITH ONLY PEAKS
    nmx=npks;
    tmp=0.*im;
    for i=1:nmx
        tmp(mx(i,1),mx(i,2))=im(mx(i,1),mx(i,2));
    end
    %LOOK IN NEIGHBORHOOD AROUND EACH PEAK, PICK THE BRIGHTEST
    for i=1:nmx
        roi=tmp( (mx(i,1)-floor(sz/2)):(mx(i,1)+(floor(sz/2)+1)),(mx(i,2)-floor(sz/2)):(mx(i,2)+(floor(sz/2)+1))) ;
        [mv,indi]=max(roi);
        [mv,indj]=max(mv);
        tmp( (mx(i,1)-floor(sz/2)):(mx(i,1)+(floor(sz/2)+1)),(mx(i,2)-floor(sz/2)):(mx(i,2)+(floor(sz/2)+1)))=0;
        tmp(mx(i,1)-floor(sz/2)+indi(indj)-1,mx(i,2)-floor(sz/2)+indj-1)=mv;
    end
    ind=find(tmp>0);
    mx=[mod(ind,nr),floor(ind/nr)+1];
end

if (size(mx)==[0,0])
    out      = [];
    out_ends = [];
    goodSpot(rr,2) = 0;  % Added 9/10/2017 gh
    fprintf('goodSpot(rr,2) = 0');
    
else
%    out(:,1) = mx(:,1); % changed 9/10/2017 gh 
%    out(:,2) = mx(:,2);
out(:,2) = mx(:,1); % restored original code June 4 2018 gh. 
out(:,1) = mx(:,2);
end
