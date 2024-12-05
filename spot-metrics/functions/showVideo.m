function [] = showVideo(stack, framerate, filename, motion)
%SHOWVIDEO Displays and image stack cell array as a movie

    % stack     cell array containing images in the stack
    % framreate rate at which frames should be displayed
    % filename  file to which movie will be written instead of shown
    % motion:   if set, shows diff to first frame for every frame
    
    if nargin < 4
        motion = false;
    end
    
    % create movie cell array
    num_frames = size(stack, 1);
    video = cell(num_frames,1);
    
    % set default framerate if none passed
    if(nargin < 2)
        framerate = 10;
    end
    
    % prepare video writer if passed
    if nargin >= 3
        outputVideo = VideoWriter(filename, 'Uncompressed AVI');
        outputVideo.FrameRate = framerate;
        open(outputVideo);
    end
    
    % populate movie with stack frames
    for f = 1 : num_frames
        
        % get current frame from stack
        cframe = stack{f};
        
        % if the content is grayscale
        if size(cframe,3) == 1
        % adjust frame for contrast
        cframe = imadjust(cframe);
        % do necessary conversions
        cframe = im2double(cframe);
        
            if f == 1
                firstframe = cframe;
            end
            
            if motion
                cframe = imfuse(cframe,firstframe);
            end
                
        
        end
        
        % write frame to video
        if nargin >= 3
            writeVideo(outputVideo,cframe);
        else
            if ~motion
                cframe = cat(3, cframe, cframe, cframe);
            end
            video(f) = im2frame(cframe); 
        end
            
    end
    
    % close video file
    if nargin >= 3
        close(outputVideo);
    else
        % set up the viewer
        figure('Name','Video','NumberTitle','off');
        imshow(video(1).cdata, 'Border', 'tight');

        % display the video as a loop
        movie(video,10,framerate);
    end
    
end
