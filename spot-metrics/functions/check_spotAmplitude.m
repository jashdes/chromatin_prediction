function [ goodSpot ] = check_spotAmplitude( im,th_percent,rr,Imax_im_filtered,goodSpot)
%inputs
    % im = ROIdata
    % th_percent = threshold in % of Imax_im_filtered
    % rr = spot number;
    % Imax_im_filtered = mzx(max(im_filtered(:,:,rr,kk=1)
    % goodSpot = 49 x 5 array 
%outputs
    % goodSpot modified
    
    if (max(max(im)) < th_percent*Imax_im_filtered)
        goodSpot(rr,2)= 0;
    else goodSpot(rr,2)=1;
    end
end


