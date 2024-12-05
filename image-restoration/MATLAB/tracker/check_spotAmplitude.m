function [ goodSpot, Imax_rr] = check_spotAmplitude( im,th_amplitude,rr,Imax_im_filtered,goodSpot,Imax_rr)
%inputs
    % im = ROIdata for one spot
    % th_amplitude = threshold as fraction of Imax_im_filtered.
    % rr = spot number;
    % Imax_im_filtered = max(max(im_filtered(:,:,rr,kk=1). Note standard 
    %       is same spot, first frame.
    % goodSpot = 49 x 7 array 
%outputs
    % goodSpot modified in column 2 and 7.
    Imax_rr(rr) = max(max(im));
    if (max(max(im)) < (th_amplitude*Imax_im_filtered))
        goodSpot(rr,2) = 0;
        goodSpot(rr,7) = 0;
    else goodSpot(rr,2)=1;
    end
end


