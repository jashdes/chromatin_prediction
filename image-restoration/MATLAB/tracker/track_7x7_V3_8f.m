function [ out_batch, SS, goodSpot] = track_7x7_V3_8f(im_array,template,stackNameFull,stackName,stackName_noTif,stackFolder,fileNameMatlabStruct,numFrames,ff,PT1,PT2, pix_per_micron, sec_per_frame);
% fileName = track_7x7_V3_8  called on L152 in Main_Batch_7x7
% folderName = F:\Userdata\George\Matlab\PAGFP_7x7_V3.8\
% computer   = HPZ4A, now at home
% GHolzwarth 2018-2020

% function calls:
    % pkfnd_   several versions in different places. 
    % pkfnd_1_spot
    % bpass
    % centerOfMass
    % signalToNoise,
    % profileFunction,
    % displayFigs_2peaks_realData(      )

% STEP 10  select scale factors, set stack YES-NO, choose directory for images
fprintf('STEP 10, track_7x7_V3.8, set scale factors, select data\n');
        % pix_per_micron = 9.36;     % Generate_movie_2_gaussians code
        %pix_per_micron = 15.18;  % Zeiss 710 datasets 2015_09_03, _09_16,l_09_28
        % pix_per_micron = 9.36;  % generate_movie 3 x 3 Gaussians  
        % pix_per_micron = 19.63;  % Olympus FV1200 confocal (dataset 2015_07)
        % stackYesNo = 1; Is data a stack or a set of tifs?
    numExpts   = 10;    % needed for msd_out file ????? 
    numSpotsX = 7;      % number of spots along x direction
    numSpotsY = 7;      % number of spots along y direction
    numSpots   = (numSpotsX)*(numSpotsY);
    numCols_msd_out = 3;  % needed for msd_out file

% STEP 11.1 obtain metadata from image 1 and construct fileNames for output
% fprintf('STEP 11.1 obtain metadata\n');
% %   ****Start
% fds = fileDatastore('*.tif', 'ReadFcn', @importdata);
% stackNameFull = fds.Files;
% numStacks = length(stackNameFull);
%         for fff = 1 : numStacks
%           [stackFolder,stackName_noTif{fff},ext] = fileparts(stackNameFull{fff})
%           stackName{fff} = strcat(stackName_noTif{fff},'.tif');
%         end
% info=imfinfo(stackNameFull{ff});
% info=imfinfo(stackNameFull);
% numRows   = info.Height;
% numCols   = info.Width;
% [numFrames,dummy] = size(info);   
info=imfinfo([stackFolder stackName]);
numRows   = info.Height;   % 181 for 7x7 Olympus
numCols   = info.Width;    % 181   "
% parts     = regexp(stackNameFull, '\', 'split');
% DirPart   = char(parts(end-1));
fileNameExcel  = strcat(stackName_noTif,'10frames.xlsx');
fileNameMatlab = strcat(stackName_noTif,'_msd_10frames.mat');
fileNameCSV    = strcat(stackName_noTif,'_msd_DIS.csv');
fileNameMatlabStructA = strcat(stackName_noTif,'_track_S_A.mat');
fileNameMatlabStruct = strcat(stackName_noTif,'_track_S.mat');
rr_max = numSpots;
kk_max = numFrames;
        
% STEP 11.2  Set width, define ROI, set width2, define ROI2
width = 8;  % use with template-based ROI centers
numColsROI = 2*width + 1;   % 17
numRowsROI = 2*width + 1;
width2 = 4; % use with CL centered ROI2
numColsROI2 = 2*width2 + 1;         % 9
numRowsROI2 = 2*width2 + 1;         % 9

% STEP 11.3 define t_array and t_array_sec
t_array        = linspace(0, (numFrames -1),numFrames); % frameNumber
t_array_sec    = t_array*sec_per_frame;                 % frameTime

% STEP 12.1 construct pedestal for each frame 
fprintf('STEP 12.1 pedestal construction\n');
pedestal = zeros(numRows,numCols,numFrames);
im_minus_pedestal   = zeros(numRows,numCols,numFrames);%       "
SE = offsetstrel('ball',5,10);
     %SE = offsetstrel('ball',R,H) creates a nonflat, ball-shaped
     % structuring element whose radius in the X-Y plane is R and whose
     % maximum offset height is H. For improved performance,
     % offsetstrel approximates this shape by a sequence of eight nonflat
     % line-shaped structuring elements.
  
for kk = 1:numFrames
     pedestal(:,:,kk) = imopen(im_array(:,:,kk), SE);
end  % end for kk = 1;numFrames

% STEP 12.2 Smooth pedestal in time domain
sigma = 2;
[pedestalFilt] = gaussfilt(pedestal, sigma);
    
% STEP 12.3 smooth pedestal over time and fix problems at ends
for kkk = 1:(3*sigma -1)
   pedestalFilt(:,:,kkk) =   pedestalFilt(:,:,3*sigma);  % left end fix
   pedestalFilt(:,:,(end - (kkk-1)))= pedestalFilt(:,:,(end - (3*sigma-1))); % right end fix
end
for kk = 1:numFrames  % redo ypeak
   ypeak(kk) = max(max(pedestalFilt(:,:,kk)));
end
    
% STEP 12.4 Subtract pedestalFilt from im_array
Imax_kk=zeros(1,numFrames);
for kk = 1:numFrames
        im_minus_pedestal(:,:,kk) = im_array(:,:,kk) - pedestalFilt(:,:,kk);
        Imax_kk(kk)=max(max(im_minus_pedestal(:,:,kk)));
        if (im_minus_pedestal(:,:,kk)<0)
               im_minus_pedestal(:,:,kk)=0;
        end    % end if
end  % end for kk = 1;numFrames

% STEP 12.5 Display improfile for im, pedestal, and im - pedestal,
            % for kk = 1 and kk = end
       
kk = 1;      
x = [1, 181];  % beginning and end of line profile
y = [91, 91];  % exact center pixel of 181 x 181
ymax = max(max(im_array(:,:,kk))) + 100;

% figure (10)  %  *******   profile, background, profile-backgd ***  
% subplot(2,3,1)
%   linprof_im  = improfile(im_array(:,:,1),x,y);
%    plot(linprof_im);
%      xlim([1 181]);
%      ylim([0 ymax]);
%      title('10A image 1');
% 
% subplot(2,3,2)
%    linprof_pedestal = improfile(pedestal(:,:,1),x,y);
%    plot(linprof_pedestal);
%      xlim([1 181]);
%      ylim([0 ymax]);
%      title('10B bkgd ');
%     
% subplot(2,3,3)
%    linprof_sub = improfile(im_minus_pedestal(:,:,1),x,y);
%    plot(linprof_sub);
%      xlim([1 181]);
%      ylim([0 ymax]);
%      title('10C im-bkgd');
% 
%      kk = numFrames
% subplot(2,3,4)
%   linprof_im2  = improfile(im_array(:,:,kk),x,y);
%    plot(linprof_im2);
%      xlim([1 181]);
%      ylim([0 ymax]);
%      title('10D last image ');
% 
% subplot(2,3,5)
%    linprof_pedestal2 = improfile(pedestal(:,:,kk),x,y);
%    plot(linprof_pedestal2);
%      xlim([1 181]);
%      ylim([0 ymax]);
%      title('10E bkgd ');
%     
% subplot(2,3,6)
%    linprof_sub2 = improfile(im_minus_pedestal(:,:,kk),x,y);
%    plot(linprof_sub2);
%      xlim([1 181]);
%      ylim([0 ymax]);
%      title('10F im-bkgd2');
          
% STEP 14 Evaluate the maximum intensity in im_minus_pedestal stack.
Imax_im_minus_pedestal = max(max(max(im_minus_pedestal))); 

% STEP 15  Apply bpass filter to im-minus-pedestal
fprintf('STEP 13 bpass\n');
for kk = 1:numFrames
    im_filtered(:,:,kk) = bpass_plusEdges(im_minus_pedestal(:,:,kk),1,10,0.05);  % *****     bpass
end  % loop over kk, L202-236

% STEP 16 set up registration
fprintf('STEP 16 registration\n');

% STEP 16.1 Scale intensities of template and image to 1.0   
    Imax_template    = max(max(template(:,:,1)));
       template_scaled = template./Imax_template;
    Imax_im_filtered = max(max(im_filtered(:,:,1))); % 2000 for BLEO 10
       im_filtered_scaled = im_filtered./Imax_im_filtered;

%  STEP 16.2 register
[optimizer, metric]  = imregconfig('multimodal');
optimizer.MaximumIterations = 500; 
optimizer.InitialRadius = 0.002;
 
[templateRegistered Rfixed] = imregister(template_scaled(:,:,1),im_filtered_scaled(:,:,1),'affine',optimizer,metric,'DisplayOptimization',0);

% STEP 16.3 Find tform and apply it to template
tform = imregtform(template_scaled(:,:,1),im_filtered_scaled(:,:,1),'affine',optimizer,metric,'DisplayOptimization',0);
       tform.T
       tform_array = maketform('affine',tform.T)
       
% STEP 16.35 Hard-code the original centers for 7x7 template in 181 x 181 image
fprintf('STEP 15.1 hard-code 7x7 centers.');
  %  Hard-code centers based on pattern generated by DOE in Olympus
        % using 181 x 181 ROI with center spot at (91,91) and spacing =
        % 16.67. Check! this doesn't compute correctly
C=zeros(49,2);  %col2=1, x; col2 = 2, y

C(:,1) = [44 44 44 44 44 44 44 60 60 60 60 60 60 60 75 75 75 75 75 75 75 91 91 91 91 91 91 91 107 107 107 107 107 107 107 122 122 122 122 122 122 122 138 138 138 138 138 138 138];
C(:,2) = [44 60 75 91 107 122 138 44 60 75 91 107 122 138 44 60 75 91 107 122 138 44 60 75 91 107 122 138 44 60 75 91 107 122 138 44 60 75 91 107 122 138 44 60 75 91 107 122 138];
u = C(:,1)';    % u is 1 x 49 array  Note TRANSFORM. C(:,1) gives columns (x)
v = C(:,2)';    % v     "            Note TRANSFORM. C(:,2) gives rows    (y)
% transform the centers
[xm,ym] = tformfwd(tform_array,u,v)  % xm is 1 x 49 array of x values of centers
                                     % use xm in place of C(:,1)
                                     % use ym in place of C(:,2)
xmr = round(xm + 1);  % +1 is fudge factor
ymr = round(ym + 1);
    
% STEP 16.4 construct 49 masks to assign a number to each spot
fprintf('STEP 15.2 Construct ROI-s and create goodSpot\n');
rr = linspace(1,49,49);
M     = zeros(numRows,numCols,49);
mask0 = zeros(numRows,numCols,49);
mask  = uint8(zeros(numRows,numCols,49));
imMask = zeros(numRows,numCols,49);
templateRegMasked = zeros(numRows,numCols,49);  % 
im_filtMasked     = zeros(numRows,numCols,49);  % new June 2017

for rr = 1:49 % end at L127
   M((ymr(1,rr)-width):(ymr(1,rr) + width),(xmr(1,rr)-width):(xmr(1,rr)+width),rr) = 1; %changed 2018_06_05
        %  M((xmr(1,rr)-width):(xmr(1,rr) + width),(ymr(1,rr)-width):(ymr(1,rr)+width),rr) = 1;
   % M is 0 except in ROI
     %   M((C(rr,2)-width):(C(rr,2) + width),(C(rr,1)-width):(C(rr,1)+width),rr) = 1;
    imMask(:,:,rr) = mat2gray(M(:,:,rr));
    mask(:,:,rr)   = uint8(M(:,:,rr));  % I = mat2gray(A) converts the matrix A to the intensity image I
    mask_double(:,:,rr) = double(M(:,:,rr));      
    im_filtMasked(:,:,rr)     = im_filtered(:,:,1).*mask_double(:,:,rr);  % 181 x 181 x numFrames
  
end  % end of rr loop
fprintf('finished marking ROIs of im_filt');

% STEP 16.5 Set up goodSpot matrix
fprintf('STEP 15.3 Set up goodSpot matrix\n');
goodSpot = ones(49,7); %  col6 and 7 added 2020_04
    %col1     2           3          4       5        6        7
    % rr  Amplitude   max(resnorm)  NaN   position   dD/D   decision
    %                                       in ROI 
    % goodAmplitude = 1 means bright spot; goodAmplitude = 0 means weak spot
rr_array = linspace(1,49,49);
goodSpot(:,1) = rr_array';

% STEP 16.6 Find center of each spot using pkfnd_7x7_original
fprintf('STEP 16 find center of each spot with pkfnd_7x7_original\n');
out_temp_reg = zeros(numSpots,2);  % consider cleanup. 
Imax_rr= zeros(1,numSpots);
for rr = 1:numSpots % short rr loop ends L273
    th_amplitude = .05;  % fraction of Imax, threshold
    % goodSpot is 49 x 7.
    [goodSpot,Imax_rr] = check_spotAmplitude( im_filtMasked(:,:,rr),th_amplitude,rr,Imax_im_filtered,goodSpot,Imax_rr);        %[goodSpot(ff)] = check_spotAmplitude( im_filtMasked(:,:,rr),th_amplitude,rr,Imax_im_filtered,goodS(ff)pot);
    Imax = max(max(im_filtMasked(:,:,rr)));
    if(goodSpot(rr,2)==1)  % This function call seems to duplicate earlier stuff
        [out_no_rr goodSpot2] = pkfnd_7x7_original(im_filtMasked(:,:,rr),0.9*Imax,3,rr,goodSpot);
        % goodSpot2 should be a dead end. 
    else
        out_no_rr(1) = NaN;
        out_no_rr(2) = NaN;
    end
       
    out_temp_reg(rr,1)=out_no_rr(1,1);  % "out_temp_reg" replaces "out"
    out_temp_reg(rr,2)=out_no_rr(1,2);
end  % end of short rr loop 257:273
fprintf('L294, goodSpot next');  goodSpot(1:49, 1:7)
pause (PT1) % goodSpot correct col2 col7

% STEP 16.7 Construct square box centered on "out_temp_reg" for each rr
ROI_template = zeros(49,4);  % rows = spot#, col1 = ii_initial; col2 = ii_final;
                              % col3 = jj_initial; col4 = jj-final;
Imax = max(max(templateRegistered(:,:,1))); 

% STEP 16.75 define ROI around each spot in template
ROI_template = zeros(numSpots,4);  % rows = spot#, col1 = ii_initial; col2 = ii_final;
                    %               col3 = jj_initial; col4 = jj-final;
kk = 1;  % use first image to define all ROIs in next loop
 
for rr = 1:rr_max % rr is ROINumber L243-339
    %     ROI_template(rr,1) = round(xmr(1,rr)); % x_i along Columns
    %     ROI_template(rr,2) = round(xmr(1,rr)+2*width); % x_f along Columns 
    %     ROI_template(rr,3) = round(ymr(1,rr)); % y_i down ROWS
    %     ROI_template(rr,4) = round(ymr(1,rr)+2*width); % y_f down ROWS 
ROI_template(rr,1) = round(xmr(1,rr)-width); % x_i along Columns
ROI_template(rr,2) = round(xmr(1,rr)+width); % x_f along Columns 
ROI_template(rr,3) = round(ymr(1,rr)-width); % y_i down ROWS
ROI_template(rr,4) = round(ymr(1,rr)+width); % y_f down ROWS 
end  % end of 'for rr = 1:rr_max' L 277-402

% STEP 16.8 construct ROI_im_filtdata
ROIdata_template =  zeros(numRowsROI,  numColsROI,  49, numFrames); % 2*width + 1
for rr = 1:numSpots
    if(goodSpot(rr,2)==1)
        row_begin = ROI_template(rr,3)  % y_i
        row_end   = ROI_template(rr,4)  % y_f
        col_begin = ROI_template(rr,1)  % x_i 
        col_end   = ROI_template(rr,2)  % x_f
        for kk = 1:numFrames
            ROIdata_template(:,:,rr,kk)= im_filtered(row_begin:row_end,col_begin:col_end,kk);
        end
    else 
        % ROIdata_template(:,:,rr,kk)= im_filtered(row_begin:row_end,col_begin:col_end,kk).*0;
    end
end  % end of rr for loop


% % Fig 2 start
% 
%     if (ff == 1)
%         fprintf('STEP 2.5 read in stack{1} and show it\n');
%             for kk=1:numFrames
%                 im_array(:,:,kk)=double(imread(stackNameFull,'Index',kk));
%                 % im_array(:,:,kk)=double(imread(stackFolder stackName{1}],'Index',kk));
%             end
%         Imax = max(max(max(im_array)));
%         hfig2 =figure (2) %********************  Fig 2 ******************* 
%         axes2 = axes('Parent',hfig2); hold(axes2, 'on');
%         xlim(axes2,[0.5 181.5]);
%         ylim(axes2,[0.5 181.5]);
%         box(axes2,'on');
%         axis(axes2,'ij');
%         set(axes2,'DataAspectRatio',[1 1 1],'Layer','top',...
%             'TickDir','out','XTick',[20 40 60 80 100 120 140 160 180],...
%             'YTick',[20 40 60 80 100 120 140 160 180]);
%        for kkk=1:numFrames
%             imm(:,:) = im_array(:,:,kkk);
%             imshow(imm,[]); % [Imax/2 Imax]
%             title(['Fig 2. ff = ',num2str(ff),' kk = ',num2str(kkk)]);
%             pause(PT1)
%        end
% fprintf('stop');

% % *****    

% Fig2 end 
    % ROIdata = zeros(17,17,numFrames);
    % for kk=1:numFrames
    %   ROIdata(:,:,kk)=double(ROIdata_template(:,:,49,kk));
    % end
    % fprintf('track L310\n')
% pause (PT2);
    % hfig70 = figure (70)
    % imshow(ROIdata_template(:,:,1,1),[]);
    %     title('Fig70 ROIdata_template(:,:,1,1)');
    %     pause ();

% STEP 17 Set up main loop to fit 2D Gaussian to each peak.
fprintf('STEP 18.0 set up main loop fitting 2D Gaussian\n');
track_out = zeros(numFrames*numSpots,13); % 200 x 49 = 9800
    % cols:  col1 = spot(rr); col2 = frame(kk); col3 = x; col4 = y; col5 = t
    % col 11=flag for strong/weak spot, col12 = sigma_x, col13=sigma_y.

% declare arrays needed inside loops.           
fitparams_kk = zeros(numFrames,5);
dataSmall      = zeros(numColsROI,numRowsROI);
dataSmall2     = zeros(numColsROI2, numRowsROI2);
CLrrkk_corr    = zeros(1,2,rr_max,numFrames); % added 2017_01_27

% % ***********************************************************************
% % ***********    MAIN LOOPS over spot=rr and frame=kk      **************
% % **********    subpixel localization of each spot     ******************
% % ***********************************************************************

% STEP 18 MAIN LOOPS fitting 2D Gaussian to each spot in each frame. 
dataSmall = zeros(numColsROI,numRowsROI);
    % goodSpot = zeros(49,5); % 
        % col    1           2                3               4            5
        %       rr  goodAmplitude    max(resnorm) noNaN    decision 
        % goodAmplitude = 1 means bright spot; goodAmplitude = 0 means weak spot    
resnorm_rrkk = zeros(49,numFrames); % residual returned by lsqnonlin

    % fprintf('STEP 18.1 breakpoint before main loops\n');
% *************   *************   ***************   **************   ******
% *************   Main loops over rr and kk         **************   ******
% *************   *************   ***************   **************   ******
for rr = 1:numSpots % rr is spot number. loop ends 483 ***** rr loop starts
    fprintf(' rr = %s  \n',num2str(rr));
 if(goodSpot(rr,2)==1) %ends 481
    for kk = 1:numFrames   %  **************kk loop ends 465
      dataSmall(:,:) = ROIdata_template(:,:,rr,kk);
                     
% STEP 18.2. Locate intensity centroid for this spot, using pkfnd_1_spot,
                % to initialize parameter (A0).
nnn = 14;  % nnn = figure number
[CL, nnn, maxIntensity] = pkfnd_1_spot(dataSmall,nnn); % uses Center-of-Mass
if (maxIntensity > 0.1*Imax_im_filtered)  % duplicates earlier steps
    goodSpot(rr,2)=1;
    CLrrkk_corr(1,rr,kk) = CL(1);   % consider elimination of CLrrkk but retain CL
    CLrrkk_corr(2,rr,kk) = CL(2);
    coordinates(1) = round(CL(2));  % x (col) of center pixel. Note 1->2 switch 
    coordinates(2) = round(CL(1));  % y (row) of center pixel
else
        % fprintf('L491 track_7x7, inside else\n');
        CL(1)=NaN;
        CL(2)=NaN;
       goodSpot(rr,2) = 0;
    end             
b = dataSmall;  % was ROIdata(:,:,kk,rr); changed to dataSmall2 Nov 1 2016 ; changed to dataSmall 4/2017
f=b; % skips bpass step. Image filtered earlier

% % ****************************************************************
% STEP 18.3 Fit 2D Gaussian to filtered peak    ***********
% % ****************************************************************dataSmall
% STEP 18.31  set up grid and ROI (patch2)on which Gaussian is evaluated.
        aa = width;  % width of ROI for Gaussian fit, 7x7 code
       [x y]=meshgrid(1:2*aa+1,1:2*aa+1);  
       f_scaled = f./(max(max(f)));
       patch2 = f_scaled; % should get rid of the "patch2" terminology. Not helpful

% STEP 18.4 BKGD SUBTRACTION, flat
% compute bkgd from 4 edges of patch2 and subtract
        meanbackground=mean([mean(patch2(1,:)) mean(patch2(2*width2+1,:)) mean(patch2(:,1)) mean(patch2(:,2*width2+1))]);
            %meanbackground=mean([mean(patch2(1,:)) mean(patch2(2*aa+1,:)) mean(patch2(:,1)) mean(patch2(:,2*aa+1))]);
        patch2=patch2-meanbackground;
%                
% define the anonymous 2D Gaussian fit function MINUS experimental data
        myfun = @(A) A(1)*exp(-(((x(:)-A(2)).^2/(2*A(3)^2))+((y(:)-A(4)).^2/(2*A(5)^2))))-patch2(:);
             % A(1) = amplitude.   patch2 is intensity of filtered image -bkgd
             % A(2) = x-center
             % A(3) = sigma x
             % A(4) = y-center
             % A(5) = sigma y
 % define the 2D Gaussian fit function alone.
    % need x(:) and y(:)
 
%      im_fit(:,:)= A(1)*exp(-(((x(:)-A(2)).^2/(2*A(3)^2))+((y(:)-A(4)).^2/(2*A(5)^2))));
            
% STEP 18.5 % find spot in frame 1
% ***   ***************   **************   **********
% ***        FIT Gaussian to spot rr in frame 1     **********
% ***   ***************   **************   **********
 if ( (kk==1)&&(goodSpot(rr,2) ==1))
           A0=[max(max(patch2));CL(2);2;CL(1);2];
           options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
             'OptimalityTolerance',1e-7,'Display','off');
           lb = []; 
           ub = []; 

    % CALL lsqnonlin **************    ******************   *************       
       [fitparams, resnorm] = lsqnonlin(myfun,A0,lb,ub,options);
       fitparams_kk(kk,:)   = fitparams(:)'; % NOTE transpose
       resnorm_rrkk(rr,kk)  = resnorm;
                        
  end  % L531 ends "if kk==1" L695
  if((kk==1)&&(goodSpot(rr,2)==0)) % spot is too dim
       goodSpot(rr,7)= 0;  % reject spot.          
  end

% STEP 18.5 find spot in remaining frames
  if ((kk>=2)&&(goodSpot(rr,2)==1))
   options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
             'OptimalityTolerance',1e-7,'Display','off');
   % CALL lsqnonlin **************    ******************   *************       
   [fitparams,resnorm] = lsqnonlin(myfun,A0,[],[],options);
     fitparams_kk(kk,:) = fitparams(:)'; % NOTE transpose
     resnorm_rrkk(rr,kk) = resnorm;
  end
  if ((kk>=2)&&(goodSpot(rr,2)==0))
     fitparams(:)=NaN;
     fitparams_kk(kk,:) = fitparams(:)';
  end
        
% STEP 18.6 Show contour plot of peak and center found by fitting.  
% fitted x vs y for one peak superposed on color image patch2(rr,kk). 
% imshow(f_scaled,[]); hold on;
%     set(gcf,'units','normalized','outerposition',[0 0 1 1 ]);
%     colormap('jet');
%     scatter(fitparams(2),fitparams(4),'ob','filled'); hold on;
%     scatter(coordinates(:,1),coordinates(:,2),'xr');
%     title(['Fig21 rr= ',num2str(rr),' kk= ',num2str(kk)]);
%         hold off;
%         pause(PT1)
% STEP 18.9 Start track_7x7 column header definitions
track_out_column_headers = cell(2,11);
track_out_column_headers{1,1}='col1';
track_out_column_headers{1,2}='col2';
track_out_column_headers{1,3}='col3';
track_out_column_headers{1,4}='col4';
track_out_column_headers{1,5}='col5';
track_out_column_headers{1,6}='col6';
track_out_column_headers{1,7}='col7';
track_out_column_headers{1,8}='col8';
track_out_column_headers{1,9}='col9';
track_out_column_headers{1,10}='col10';
track_out_column_headers{1,11}='col11';

track_out_column_headers{2,1}='rr bead';
track_out_column_headers{2,2}='kk=frame';
track_out_column_headers{2,3}='x_pix_ROI';
track_out_column_headers{2,4}='y_pix_ROI';
track_out_column_headers{2,5}='x_pix_181x181';
track_out_column_headers{2,6}='y_pix_181x181';
track_out_column_headers{2,7}='sigma_x pix';
track_out_column_headers{2,8}='sigma_y pix';
track_out_column_headers{2,9}='Io, pix';
track_out_column_headers{2,10}='t, sec';
track_out_column_headers{2,11}='goodSpot';

% STEP 19. Define track_out, the output array from track_7x7. Catenation. 
track_out(kk+(rr-1)*numFrames,1) = rr;
track_out(kk+(rr-1)*numFrames,2) = kk;
track_out(kk+(rr-1)*numFrames,3) = fitparams(2);  % x in ROI coordinates
track_out(kk+(rr-1)*numFrames,4) = fitparams(4);  % y in ROI coordinates
track_out(kk+(rr-1)*numFrames,5) = fitparams(2) + xmr(1,rr); % x in im space  2017_04_14
track_out(kk+(rr-1)*numFrames,6) = fitparams(4) + ymr(1,rr); % y in im space  2017_04_14 
track_out(kk+(rr-1)*numFrames,7) = fitparams(3);  % sigma_x in ROI coordinates
track_out(kk+(rr-1)*numFrames,8) = fitparams(5);  % sigma_y in ROI coordinates
track_out(kk+(rr-1)*numFrames,9) = fitparams(1);  % Io = peak intensity of Gaussian
track_out(kk+(rr-1)*numFrames,10) = t_array_sec(kk);  % time (s)
track_out(kk+(rr-1)*numFrames,11) = goodSpot(rr,5);  % goodSpot
   
fitparams_TF(kk,:)   = isnan(fitparams_kk(kk,:));
fitparams_no_NaN(kk) = all(fitparams_TF(kk,:)<1);
  % 1= no NaN(good); 0 = some fitparams = NaN (bad)
    end; % end of kk=1:numFrames L328:468

% STEP 19.1 goodSpot column 3 (was 4)  L466
 fitparams_no_NaN_all_kk = all(fitparams_no_NaN(:)>0);
     % 1= no NaN(good); 0 = some fitparams = NaN (bad)
   if (fitparams_no_NaN_all_kk == 1) % no NaN, all kk, this rr is good
      goodSpot(rr,3) = 1;
      goodSpot(rr,7) = goodSpot(rr,2) & goodSpot(rr,3);
      track_out(kk+(rr-1)*numFrames,11) = 1;
   end
   if (fitparams_no_NaN_all_kk == 0) % NaN present at some kk, this rr is bad
      goodSpot(rr,3) = 0; % was 4
      goodSpot(rr,7) = 0;
      track_out(kk+(rr-1)*numFrames,11) = 0; 
   end % end of if
 else
   % goodSpot(rr,7) = 1;  % no NaN. Spot ok.
 end % end of check for NaN L325:480
end; % end of BigLoop, rr = 1:numSpots L379-481 ********** 

% Calculate goodSpot(rr,4) based on maximum allowable value of resNorm_rrkk
maxresnorm = max(resnorm_rrkk,[],2); thresh_maxresnorm = 25; % 2020_04
if (maxresnorm >= thresh_maxresnorm)
    goodSpot(:,4) = 0; % **************************    set goodSpot(rr,4)
else goodSpot(:,4)= 1;
end
% STEP 19.2 add fields to structure SS.  
SS.field_stackName_noTif = stackName_noTif; % done in Main
SS.field_track_out_column_headers = track_out_column_headers;
SS.field_folderName     = stackFolder;
SS.field_track_out      = track_out
SS.field_t_array_sec    = t_array_sec;
SS.field_pix_per_micron = pix_per_micron;
SS.field_numRows        = numRows;
SS.field_numCols        = numCols;
SS.field_numFrames      = numFrames;
SS.field_numSpots       = numSpots;
SS.field_goodSpot       = goodSpot;
SS.field_xmr            = xmr; % added 2020_04_17
SS.field_ymr            = ymr; %   "
SS.field_ROI_template   = ROI_template;  % added 2020_06_08
% % % STEP 5  Create structure out_batch to save results
out_batch.templateName   = 'templateName';   
