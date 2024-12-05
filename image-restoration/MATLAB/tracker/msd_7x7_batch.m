% msd_7x7_batch.m
    % history: based on msd_3x3_batch.m 
    % requires SS in workspace, from prior execution of track_7x7_func.
    % gholzwarth plus many others previously. All errors blame gh.
    % started October 2016. 
    % msd is a script, so vars left in the workspace by
    %  MAIN_Batch_7x7 and Track_7x7 are available to this file
% fprintf('STEP 10 MSD for each peak\n');

pause on;
PT2 = 2;    % pause time #2, long

% STEPS are numbered 30:39 for msd_7x7_batch
% STEP 30 load tracking results computed in track_7x7_V3.8

fileName   = stackName_noTif;
folderName = stackFolder;
numRows    = SS.field_numRows;
numCols    = SS.field_numCols;
numFrames  = SS.field_numFrames;
numSpots   = SS.field_numSpots;
goodSpot   = SS.field_goodSpot;  % should be 49 x 7
xmr        = SS.field_xmr;
ymr        = SS.field_ymr;
ROI_template = SS.field_ROI_template;  % 49 x 4. See L 298 in tracking_7x7

% fileName   = S(ff).field_stackName_noTif;
% folderName = S(ff).field_folderName;
% numRows    = S(ff).field_numRows;
% numCols    = S(ff).field_numCols;
% numFrames  = S(ff).field_numFrames;
% numSpots   = S(ff).field_numSpots;
% goodSpot   = S(ff).field_goodSpot;  % should be 49 x 7
% xmr        = S(ff).field_xmr;
% ymr        = S(ff).field_ymr;
pix_per_micron = SS.field_pix_per_micron;
micron_per_pix = 1/pix_per_micron;
track_out  = SS.field_track_out;  % Note that SS.field_track_out will be updated in msd.
    % specifically cols 5 and 6 will be corrected for drift
    % columns in track_out are:
%    1       2        3     4     5     6   7  8  9  10  11   
%  spotID  frameID  xROI  yROI    x     y                goodSpot 
%  rr      kk       pix   pix     pix   pix  
%                                 181   181
% STEP 30.1 pull in x,y data (in pixels) within the 181 x 181 image space
xpixROI = zeros(numSpots,numFrames);  % 49 x 200.includes good and bad spots  
ypixROI = zeros(numSpots,numFrames);
xpix181 = zeros(numSpots,numFrames);  % 49 x 200.includes good and bad spots  
ypix181 = zeros(numSpots,numFrames);

% STEP    construct xpixROI, ypixROI, xpix181, ypix 181; by decatenation.
for rr  = 1:numSpots
   if(goodSpot(rr,7)==1)
      for kk = 1:numFrames
        xpixROI(rr,kk) = track_out(kk+(rr-1)*numFrames,3); % x in ROI space  2020_05_01
        ypixROI(rr,kk) = track_out(kk+(rr-1)*numFrames,4); % y in ROI space  2020_05_01   
        xpix181(rr,kk) = track_out(kk+(rr-1)*numFrames,5); % x in 181 space  2020_05_01
        ypix181(rr,kk) = track_out(kk+(rr-1)*numFrames,6); % y in 181 space  2020_05_01 
      end
   else
      xpixROI(rr,kk) = NaN;
      ypixROI(rr,kk) = NaN;
      ypix181(rr,kk) = NaN;
      xpix181(rr,kk) = NaN;
   end
end

numSpotsX = 7;      % number of spots along x direction
numSpotsY = 7;      % number of spots along y direction
numSpots   = (numSpotsX)*(numSpotsY);
numCols_msd_out = 3;  % needed for msd_out file

% STEP 30.2 use goodSpot col5 to reject spots which get
%           close to edge of ROI. Logic corrected June 2020 gh
for rr = 1:numSpots
    
    for kk = 1:numFrames
        if(goodSpot(rr,5) == 1);
            if(((xpixROI(rr,kk)>3) && (xpixROI(rr,kk)<15))&&((ypixROI(rr,kk)>3) && (ypixROI(rr,kk)<15)))
                goodSpot(rr,5) = 1;
            else
                goodSpot(rr,5) = 0;
            end;
        end
    end  % end for kk...
    if(goodSpot(rr,5) == 0)
        goodSpot(rr,7)=0;
    end
end  % end for rr...

t_array = linspace(0,numFrames-1,numFrames);  % 1x200 double

% STEP 33.5 Fig30-39 x-cor vs y-cor
figNum = 29+ff;
textbox30{1}= stackName_noTif;
hfig30 = figure (figNum);   % ******************  Fig 23
axes30 = axes('Parent',hfig30);
LL = length(ROI_template(1,1):ROI_template(1,2)); % length of ROI edge
arrayLL  = linspace(1,LL,LL);
for rr = 1:numSpots
    if (goodSpot(rr,7)==1)
        plot(xpix181(rr,:),ypix181(rr,:),'-k'); hold on;
    end
end
set(gca,'Ydir','reverse');
    axis equal;  %added 3/1/2020 gh
% show borders of ROIs. Added 2020_06_10  *****************************
for rr = 1:numSpots
    col_begin = ROI_template(rr,1)+8; % width added 6/12
    col_end   = ROI_template(rr,2)+8;
    row_begin = ROI_template(rr,3)+8;
    row_end   = ROI_template(rr,4)+8;
    
    % ROItop
        x = linspace(col_begin,col_end,16)
        y = linspace(row_begin, row_begin,16)
    plot(x,y,'-r'); hold on; 
         
% % ROIbottom
    x = linspace(col_begin,col_end,16)
    y = linspace(row_end, row_end,16)
    plot(x,y,'-g'); hold on; 

% ROIleft
    x = linspace(col_begin,col_begin,16);
    y = linspace(row_begin, row_end,16);
    plot(x,y,'-b'); hold on; 

% ROIright
    x = linspace(col_end,col_end,16);
    y = linspace(row_begin, row_end,16);
    plot(x, y,'-k'); hold on;
end  % end of rr loop
    grid on;
    title(['Fig',num2str(figNum),' 23 xpix181 vs ypix181']);
    xlabel('xpix181 (pixel)');
    ylabel('ypix181 (pixel)');
    
    xlim(axes30,[0, 181]);
    ylim(axes30,[0, 181]);
    box(axes30, 'on');
    axis(axes30,'ij');
    annotation(hfig30,'textbox',[0.225 0.77 0.77 0.143],...
   'String',textbox30,'FitBoxToText','off','EdgeColor','none');
pause (PT2)

% Save Fig 30 Tracks in boxes
% construct fileNameTracksInBoxes unique to this cell, e.g. BLM__2
    fileNameTracksInBoxes  = strcat(stackName_noTif,'_TracksInBoxes.pdf');
    full_file_name = fullfile(folderNameOutputFigures, fileNameTracksInBoxes);
    saveas(gcf,full_file_name);
    saveas(gcf,full_file_name, 'pdf');

% STEP 33.6 Correct for rotation of the nucleus.  Center found but not used
%[Xc, Yc, theta_mean, theta_mean_std] = center_of_rotation_func(x_cor, y_cor);

% STEP 33.7 Convert x-cor and y_cor from pixels to microns
x_um = xpix181.*micron_per_pix;  % 49 x 200 double
y_um = ypix181.*micron_per_pix;  % 49    "
        
% STEP 34 calculate msd for each rr
numDeltaT = floor(numFrames/4);  % defines both numTau and largest tau. numDeltaT = tau
tau_array     = t_array(2:(numDeltaT+1));
t_array_sec   = SS.field_t_array_sec;  % SS brought in with "load" at L23 
tau_array_sec = t_array_sec(2:(numDeltaT+1));

msd_um = zeros(numSpots,numDeltaT,3); %# We'll store [mean, stdev, n]
 slope_5pts = zeros(3,49);    % row 1 = rr, row2 = slope, row3 = std_slope
 slope_all_pts = zeros(3,49); %    "          "              "
 yfit          = zeros(numSpots,20); % fit to first 20 values of tau.49 x20

% Main loop for MSD  *********************
for rr = 1:numSpots  % ends at L232
    if (goodSpot(rr,7) == 1)% ends L209 
        for dt = 1:numDeltaT 

           % clear delta_x_um; clear delta_y_um; clear squaredDisplacement_um; % added 2015_10_15 gh
            delta_x_um = x_um(rr,1+dt:end) - x_um(rr,1:end-dt);
            delta_y_um = y_um(rr,1+dt:end) - y_um(rr,1:end-dt);
            % msd calculated next line. Key step   K=L172
            squaredDisplacement_um= delta_x_um.^2 + delta_y_um.^2;
            msd_um(rr,dt,1) = mean(squaredDisplacement_um);   %# average
            % Trap cases where msd = NaN
            msd_um(rr,dt,2) = std(squaredDisplacement_um);    %# std
            msd_um(rr,dt,3) = length(squaredDisplacement_um); %# n = number of readings in msd

            clear delta_x_um; clear delta_y_um; clear squaredDisplacement_um; % added 2015_10_15 gh
            %clear SS;  Clear after save L 450
        end  % end of loop over dt.
    % THIS COMPLETES calculation of msd for a particular rr 
              
    % STEP 35  fit quadratic to msd, using first 5 values of tau
        S_20 = polyfitn(tau_array_sec(1:20),msd_um(rr,1:20,1),2); 
        a1_fit = S_20.Coefficients(1);
        a2_fit = S_20.Coefficients(2);
        a3_fit = S_20.Coefficients(3);
        a1_std = S_20.ParameterStd(1);
        a2_std = S_20.ParameterStd(2);
        a3_std = S_20.ParameterStd(3);
        slopePoly_20(1,rr) = rr;
        slopePoly20(2,rr) = a2_fit;
        slopePoly20(3,rr) = a2_std;
       y1 = a1_fit.*(tau_array_sec(1:20).^2);
       y2 = a2_fit.*tau_array_sec(1:20);
       y3 = a3_fit*ones(1,20);
       yfit_poly(rr,:)= y1 + y2 + y3;  

       % evaluate D and Standard Error for polynomial fit:
       slopePoly20(2,rr) = a2_fit;
       slopePoly20(3,rr) = a2_std;
       Dpoly20(rr)    = slopePoly20(2,rr)/4;  % micrometers^2/s
       SE_Dpoly20(rr) = slopePoly20(3,rr)/4;  % micrometers^2/s
       % % end of if section controlled by goodSpot(6)==1
   else
       Dpoly20(rr)    = NaN;
       SE_Dpoly20(rr) = NaN;
   end % end of if started L174   
end      % L217 end of loop over rr. Loop started L182   *****************

% STEP 35.5  use dD/D to set goodSpot 6. 2020_05_02
for rr = 1:numSpots
    if (goodSpot(rr,7) == 1)  % ends L65
        dD_over_D(rr) =  SE_Dpoly20(rr)/Dpoly20(rr);
        if (dD_over_D > 0.5) % threshold 0.5 ****************************
            goodSpot(rr,6) = 0;
            goodSpot(rr,7) = 0; 
        end  % end inner if
    end % end outer if
end  % end for rr... 

% STEP 36. Set up figure 16 (CirclePlot)showing circles with radius
    % proportional to D.  x,y in pixels
data = zeros(49,5); % 49 = spotNum, col1= spotNum; col2=x,col3=y,col4=D,col5=MarkerSize
data(1:49,1)=linspace(1,49,49); % fills col1 = rr;
for rr = 1:49
    data(rr,2)=  xpix181(rr,1);
    data(rr,3)=  ypix181(rr,1);
        % data(rr,2)=  x_cor(rr,1);
        % data(rr,3)=  y_cor(rr,1);
    data(rr,4)=  Dpoly20(rr)*1E06; % UNITS: nm^2/s
end
bins =   [ 0.0  .05  .085  .125 .175 .225  .275  .325 .375 .425 0.475 .525 .575 .625 .675 .725 ]*1E03; %units nm^2/s
bin_mean=[  0.033 0.067 0.10  0.15 0.20  0.25  0.30  0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 ]*1E03; %units nm^2/s
MarkerSize = [ 3    5     7    9   11   13   16   18   20   23   25   28   31   34   36 ];
colormap('jet');
ccc = jet(15); %divide intensity scale into 15 levels
MarkerStr = cell(1,49); % array of 49 strings 
MarkerCol = zeros(49,3); % 49x3 double for MarkerColor (r,g,b).

% Assign MarkerStr, MarkerSize, MarkerCol to each rr according to its D
for rr = 1:49   % include only spots with goodSpot = 1
    if data(rr,4)<0;  % negative D.....bad or noisy data?
       k=1;
       MarkerStr{rr} = 'o';
       data(rr,5)=1;  % default to avoid k=[] if D<0
    elseif (( data(rr,4)>0) & (data(rr,4)<=bins(end)))
       k = find(data(rr,4)>bins(1:end-1) & data(rr,4)<=bins(2:end), 1);
            %  elseif (( data(rr,4)>0) & (3*data(rr,4)<=bins(end)))
            %    k = find(3*data(rr,4)>bins(1:end-1) & 3*data(rr,4)<=bins(2:end), 1);
       MarkerStr{rr} = 'o';
       data(rr,5)=MarkerSize(k);    % 49 x 5
       MarkerCol(rr,:) = ccc(k,:);  % 15 x 3
    else 
       MarkerStr{rr} = 'x';
       data(rr,5) = 24; % oversize D . . . bad data? dropped transpose 5/13
    end
end

% STEP 36.1. Plot CircleFig with x,y in microns and quantitative colorscale.
for rr = 1:49
            data(rr,2)=  x_um(rr,1);  % x in first image
            data(rr,3)=  y_um(rr,1);  % y in first image
            data(rr,4)=  Dpoly20(rr); % Diffusion coeff using 20 time points
end
figNum = 15+ff;
textbox16{1}= stackName_noTif;        
hfig16 = figure (figNum);  %  Circle plot xy in microns ****** Fig 16,17,
% axes #1  left, bottom, width, height
   axes1 = axes('Parent', hfig16,'Position',[0.13  0.11 0.57 0.85]);
   hold(axes1,'on');
   set(gca, 'Ydir', 'reverse');
   xlim(axes1,[0, 20]); % microns
   ylim(axes1,[0, 20]); % microns
   daspect([1 1 1]);
box(axes1, 'on');
axis(axes1,'ij');
for rr = 1:49
    if (goodSpot(rr,7)==1)  % was 6 May 15 2020
        plot(data(rr,2),data(rr,3),'Marker',MarkerStr{rr},'MarkerSize',data(rr,5),'LineStyle','none',...
            'MarkerFaceColor', MarkerCol(rr,:),'MarkerEdgeColor','k');
    end
end;
xlabel('\mum');
ylabel('\mum');
annotation(hfig16,'textbox',[0.225 0.77 0.77 0.143],...
   'String',textbox16,'FitBoxToText','off','EdgeColor','none');

% axes #2   
   axes2 = axes('Parent', hfig16,'Position',[0.70  0.152 0.159 0.75]);
hold(axes2,'on');
MarkerStr2 = cell(1,15) % array of 15 strings 
MarkerCol2 = zeros(15,3); % 15x3 double for MarkerColor2 (r,g,b).
data2 = zeros(15,5);  % col1 = rr; col2=x; col3=y, col4=D, col5 = MarkerSize;
data2(:,2) = 0.9 ;               % x
data2(:,3) = linspace(2,13,15);  % y
data2(:,4) = [ 33, 67, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700 ]; % D in nm2/sec

for ss = 1:14
   fprintf('ss = '); ss;
   [row col] = find(data2(ss,4)>bins(1:end-1) & data2(ss,4)<=bins(2:end), 1);
   row
   col
   k = col;
   data2(ss,5)=MarkerSize(1,k);
   MarkerCol2(ss,:) = ccc(ss,:);  % 15 x 3
 end
 textbox15{1}= 'D (nm^{2}/s)';
 textbox33{1} = '-----33';
 textbox100{1} = '----100'; 
 textbox200{1} = '----200';
 textbox300{1} = '---300';
 textbox400{1} = '--400';
 textbox500{1} = '-500';
 textbox600{1} = '-600'; 

hold(axes2,'on');
    xlim([-0.5 3]);
    ylim([1  15]);
    set (axes2,'XTick', [],'YTick',[]);
for ss = 1:13
  plot(data2(ss,2),data2(ss,3),'Marker','o','MarkerSize',data2(ss,5),'LineStyle','none','MarkerFaceColor',MarkerCol2(ss,:),'MarkerEdgeColor','k');
  hold(axes2,'on');
end
box(axes2,'off');
set(axes2,'XTick',zeros(1,0),'XTickLabel',{},'YTick',zeros(1,0),'YTickLabel',{});
 
% ANNOTATIONS for SCALE     *****
FS = 12;  % Set FontSize globally
% 'D (nm^{2}/s)'  textbox15;
    annotation(hfig16,'textbox',[0.733 .702 0.77 0.143],...
 'String',textbox15,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
% 600 textbox600
   annotation(hfig16,'textbox',[0.80 0.609 0.95 0.143],...
 'String',textbox600,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
%  500  textbox500
   annotation(hfig16,'textbox',[0.801 0.522 0.9 0.143],...
 'String',textbox500,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
%  400 textbox400
  annotation(hfig16,'textbox',[0.790 0.436 0.77 0.143],...
  'String',textbox400,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
%  300 textbox300
 annotation(hfig16,'textbox',[0.780 0.352  0.77 0.143],...
  'String',textbox300,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
%  200 textbox200
 annotation(hfig16,'textbox',[0.771 0.268  0.77 0.143],...
  'String',textbox200,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
%  100 textbox100          
annotation(hfig16,'textbox',[0.771 0.188  0.77 0.143],...
  'String',textbox100,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
%  33 textbox33
 annotation(hfig16,'textbox',[0.772 0.103 0.77 0.143],...
  'String',textbox33,'FitBoxToText','off','EdgeColor','none','FontSize',FS);
 pause (3)
  
% save Circle_Figure (Fig16)
        % construct fileNameDfig unique to this cell, e.g. BLM__2
        % stackName_noTif
        % stackName_noTif{ff}
        % stackName_noTif(ff)
        % [pathstr2, fileName_noTif2, ext2]= fileparts([folderName fileName]);    
    fileNameCirclesFig  = strcat(stackName_noTif,'_Circle_um.pdf');
        % pathstr2_full = [pathstr2 '\CircleFigureFolder\'];
        % full_file_name = [pathstr2_full fileNameCirclesFig];
    full_file_name = fullfile(folderNameOutputFigures,fileNameCirclesFig);% 2019_12_25
    saveas(gcf,full_file_name);
    saveas(gcf, full_file_name, 'pdf');
% % % END  Fig 16:25
% 
% ***********************************************************
% STEP 12.3 compute msd_mean over 49 spots in one cell
msd_mean = mean(msd_um(:,:,1),1);  
msd_mean_std = std(msd_um(:,:,1),1);
 % ***********************
 %  % fit quadratic to msd_mean, using first 20 values of tau
 % fprintf('find slope, 2nd order polynom, 20 points\n');
     S_20 = polyfitn(tau_array_sec(1:20),msd_mean(1,1:20),2);
         %S_all = polyfitn(tau_array_sec,msd_um(rr,:,1),2);
     a1_fit = S_20.Coefficients(1);
     a2_fit = S_20.Coefficients(2);
     a3_fit = S_20.Coefficients(3);
     a1_std = S_20.ParameterStd(1);
     a2_std = S_20.ParameterStd(2);
     a3_std = S_20.ParameterStd(3);
     slopePoly_20(1,rr) = rr;
     slopePoly20(2,rr) = a2_fit;
     slopePoly20(3,rr) = a2_std;
    y1 = a1_fit.*(tau_array_sec(1:20).^2);
    y2 = a2_fit.*tau_array_sec(1:20);
    y3 = a3_fit*ones(1,20);
    yfit_poly_mean(:)= y1 + y2 + y3;
%    
% evaluate Dmean and Standard Error of Dmean from polynomial fit:
% fprintf('D from fit of polynomial to first 20 values of tau, using polyfitn\n');
slopePoly20mean(2) = a2_fit;
slopePoly20mean(3) = a2_std;
Dpoly20mean    = slopePoly20mean(2)/4;  % This is Dmean in
SE_Dpoly20mean = slopePoly20mean(3)/4;
 
% STEP 38  Add msd, D, dD, Dmean, dDmean to S
SS.field_msd_um(:,:) = msd_um(:,:,1); % col 1 = msd).
SS.field_D  =    (Dpoly20.*1e06)';  %  D in nm^2/s Note transpose
SS.field_dD = (SE_Dpoly20.*1e06)';  % dD in nm^2/s  Note transpose
SS.field_Dmean = Dpoly20mean*1e06;% Dmean in nm^2/s mean over 49 spots.
    Dmean = Dpoly20mean*1e06
SS.field_SE_Dmean = SE_Dpoly20mean.*1e06'; % std of Dmean in nm^2/s
    % SS.field_driftCorFlag = driftCorFlag;

        % %STEP 38.5 Compare diffusion distance to drift distance
        % % for at tau = t = 3s. Added 5/29 2020
        % % Diffusion distance
        %      Ldiff = (4*Dmean*3)^(0.5); %  units are nm if D is in nm^2/s
        % % Drift distance
        %      Ldrift_pix = (Drift(10,1)^2 + Drift(10,2)^2 )^(0.5); % 10th data point 
        %      Ldrift_um  =  Ldrift_pix/9.26;   % pix_per_micron = 9.26  
        %      Ldrift_nm  = Ldrift_um*1e03;     % 1000 nm/um
        %      % corresponds to 3 sec, approx
        % % Ratio     
        %      Ldiff_over_Ldrift = Ldiff/Ldrift_nm;
        %      SS.field_Ldiff_over_Ldrift = Ldiff/Ldrift_nm;
        %      % 

     
% STEP 39 update track_out cols 5,6 using drift-corrected values of x and y 
for rr = 1:numSpots
    for kk = 1:numFrames   %  **************kk loop ends 607
        track_out(kk+(rr-1)*numFrames,5) = xpix181(rr,kk); % x
        track_out(kk+(rr-1)*numFrames,6) = ypix181(rr,kk); % y 
        % track_out(kk+(rr-1)*numFrames,5) = x_cor(rr,kk) + xmr(1,rr); % x
        % track_out(kk+(rr-1)*numFrames,6) = y_cor(rr,kk) + ymr(1,rr); % y 
    end
end;

% STEP 39.1 update SS and save to folder BATCH_outputData
SS.field_track_out= track_out;
        % fileName = 'fileNameMatlabStruct'  %name of saved mat file.
        % fileName = 'fullFileNameOutputSS.mat'
save(fullfile(folderNameOutputS,fileNameMatlabStruct),'SS')
% clear SS and reset working directory 
    ffff = fields(SS);
    for kkkk=1:numel(ffff)
      SS.(ffff{kkkk})=[];
    end
% *************************************************************************
