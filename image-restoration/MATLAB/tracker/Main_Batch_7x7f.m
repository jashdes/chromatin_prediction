% Main_Batch_7x7.m    V3.8  June 7, 2020
% saved on Olin 213 HP Z4A desktop now at GH home
% folderName = F:\Userdata\George\Matlab\PAGFP_7x7_V3.8\
% based on PAGFP_7x7_V3.3PV 
% George Holzwarth 2016-2020
% 2 primary files called:
    % track_7x7_V3_8f
    % msd_7x7_batch
% 7 secondary files called:
    % bpass_plusEdges, pkfnd_7x7_original, pkfnd_modified,
    % pkfnd_1_spot, center_of_mass, center_of_rotation,
    % polyfitn
% images needed
    % tif stack(s) with PAGFP images of cells
    % tif template, synthetic data with 49 2D Gaussians spaced exactly as
        % in PAGFP data. Template makes it easier to find the 49 spots
        
function [ out_batch, SS, goodSpot] = Main_Batch_7x7f(dataFolderName,folderName, template, pix_per_micron,sec_per_frame, numFrames)
                

% STEP 0.0 Clean up workspace 
tic()
pause on;
PT1 = 0.3; % short pause time
PT2 = 5;    % long pause time

% STEP 0.1  First time setup:
% Create 5 folders on your computer to hold input and output data.
% On gholz HP Z4 A computer these folders and subfolders are: 
     % folderName = F:\Userdata\George\Matlab\PAGFP_7x7_V3.8
        % BATCH_template      = subfolder for template
        % BATCH_inputData     = subfolder for input data (tif stacks)
        % BATCH_outputFigures = subfolder for output figures
        % BATCH_outputSS      = subfolder for output structure SS.mat etc
        % BATCH_outputExcel   = subfolder for Excel spreadsheet
%
% STEP 0.2 with every new experiment
    % update template file in BATCH_template folder, if necessary 
    % put new image stacks into BATCH_inputData folder
    % modify output Excel sheetname to reflect experiment date and
    % experimental conditions.

% STEP 1.0. Get templateFileName, templateFolderName
template = imread(template);
% STEP 1.3. read and show template
figure (1) % *******************  
imshow(template);
title('Fig 1 7x7 template');
         
% STEP 2. get stackNameFull, stackName, stackName_noTif, stackFolder
dataFolderName = strcat(dataFolderName, filesep);
dataFileNameList = dir([dataFolderName '*.tif']);  % 4 x 1 cell with 6 fields
dataFileName = dataFileNameList(1).name;
stackName = dataFileName;  % ff = 1
[numStacks dummy]=size(dataFileNameList);

% stackName = dataFileNameList.name;
     % dataFileName = 'cell4.tif';
     % dataFolderName = 'F:\.......\BATCH_inputData
     % dataFileNameList = 4 x 1 struct with 6 fields.
% for ff = 1 : numStacks
%   [stackFolder,stackName_noTif{ff},ext] = fileparts(stackNameFull{ff});
%   stackName{ff} = strcat(stackName_noTif{ff},'.tif');
%   S(ff).field_stackName_noTif = stackName_noTif(ff); % S is numStacks x 1
% end

info = imfinfo([dataFolderName dataFileName]);
numRows   = info.Height;
numCols   = info.Width;
if numFrames == 0
    [numFrames,dummy] = size(info); 
end
    % numFrames  = 200; % Overide
excel_matrix         = zeros(49, 1 + 2*numStacks);  % was zeros(9, 1 + 2*numStacks);
excel_matrix(1:49,1) = linspace(1,49,49);  % rr values in first col of excel matrix
fileName_matrix      = cell(1,numStacks); % file names for excel col headings
im_array = double(zeros(numRows,numCols,numFrames));  % 181 x 181 x 200. 

% STEP 1.4 Construct output folderNames
folderNameOutputFigures = fullfile(folderName, 'BATCH_outputFigures');
mkdir(folderNameOutputFigures);
folderNameOutputS       = fullfile(folderName, 'BATCH_outputSS');
mkdir(folderNameOutputS);
folderNameOutputExcel   = fullfile(folderName, 'BATCH_outputExcel');% OK
mkdir(folderNameOutputExcel);
% 
%   ********************* *********************** **********************
%   *********************   Start Loop over ff    ********************** 
% STEP 2.1 Loop over ff, the stack index          ********************** 
fprintf('STEP 2.1 ff loop ****************\n');

for ff = 1:numStacks  % ff loop starts L110 ends 187
    try
    % **************************************************************
    % **************************************************************
    % **************************************************************
  dataFileName = dataFileNameList(ff).name;
  stackName    = dataFileName;
        stackNameFF{ff} = dataFileName; 
  [pathstr,stackName_noTif,ext] = fileparts(stackName);
        [pathstr,stackNameFF_noTif{ff},ext] = fileparts(stackNameFF{ff});
        fullFileNameOutputExcel = fullfile(folderNameOutputExcel, strcat(stackNameFF_noTif{1},'.xlsx'));
  stackFolder = dataFolderName;
  stackNameFull = fullfile(stackFolder,stackName); 
  stackNameFFFull{ff} = fullfile(stackFolder,stackNameFF{ff});
  fileNameMatlabStruct = strcat(stackName_noTif,'_SS.mat');
      if (ff == 1)
        fprintf('STEP 2.5 read in stack{1} and show it\n');
            for kk=1:numFrames
                im_array(:,:,kk)=double(imread(stackNameFull,'Index',kk));   
            end
        Imax = max(max(max(im_array)));
        hfig2 =figure (2) %*******  Fig 2 *******
        axes2 = axes('Parent',hfig2); hold(axes2, 'on');
        xlim(axes2,[0.5 181.5]);
        ylim(axes2,[0.5 181.5]);
        box(axes2,'on');
        axis(axes2,'ij');
        set(axes2,'DataAspectRatio',[1 1 1],'Layer','top',...
            'TickDir','out','XTick',[20 40 60 80 100 120 140 160 180],...
            'YTick',[20 40 60 80 100 120 140 160 180]);
       for kkk=1:10:numFrames
            imm(:,:) = im_array(:,:,kkk);
            imshow(imm,[]); % [Imax/2 Imax]
            title(['Fig 2. ff = ',num2str(ff),' kk = ',num2str(kkk)]);
            pause(PT1)
       end
    end % end if ff == 1 started L155

    % **************************************************************
    % **************************************************************
    % **************************************************************
    % STEP 2.2 ff ~=1   ************************************************** 
    if (ff~=1)
       close all; 
       for kk=1:numFrames
          im_array(:,:,kk)=double(imread([dataFolderName dataFileNameList(ff).name],'Index',kk));
       end
            display(ff);
            display(dataFileNameList(ff).name);
       figNum = 1+ff;
       %figure (figNum) %*****  Fig 3,4 ********** 
        hfig2 =figure (2) %*******  Fig 2 *******
        axes2 = axes('Parent',hfig2); hold(axes2, 'on');
        xlim(axes2,[0.5 181.5]);
        ylim(axes2,[0.5 181.5]);
        box(axes2,'on');
        axis(axes2,'ij');
       
 for kkk=1:10:numFrames
          imm(:,:) = im_array(:,:,kkk);
          imshow(imm,[]); % [Imax/2 Imax]
          title(['Fig',num2str(figNum),' ff = ',num2str(ff),' kk = ',num2str(kkk)]);
            pause (PT1)
       end
    end % end ff~=1
     
% STEP 3. track the spots
fprintf(['L155 Main kk = ',num2str(ff),'track next\n']);
[ out_batch,SS,goodSpot] = track_7x7_V3_8f(im_array,template,stackNameFull, stackName,stackName_noTif,stackFolder,fileNameMatlabStruct,numFrames,ff,PT1,PT2, pix_per_micron, sec_per_frame);

% STEP 4.  Compute msd and D for each spot, each data set.      
fprintf('STEP 4 MAIN. run msd_7x7_batch to compute msd,D for each spot\n');
    msd_7x7_batch;

% STEP 5 in Main. Finish up.  Save rotational values
% STEP 5.1 
    theta_mean = 0.0;   % XXX
    theta_mean_std = 0.0;   % XXX
    rotMatrix(ff,1) = theta_mean;
    rotMatrix(ff,2) = theta_mean_std;

% STEP 6  Put rr, D, dD, dD/D into excel_matrix   
 fileName_matrix{ff} = [fileName];
 excel_matrix(:,((ff-1)*6)+1) = linspace(1,49,49);  % rr values ff=1
 excel_matrix(:,((ff-1)*6)+2) = (Dpoly20(1,:)).*1e06;              % nm^2/s, D
 excel_matrix(:,((ff-1)*6)+3) = (SE_Dpoly20(1,:)).*1e06;           % nm^2/s; dD
 excel_matrix(:,((ff-1)*6)+4) = (SE_Dpoly20(1,:))./(Dpoly20(1,:)); % dD/D
 excel_matrix(:,((ff-1)*6)+5) = NaN;
 excel_matrix(:,((ff-1)*6)+6) = NaN;
    catch err
        warning('An error occured while processing %s\n%s\n', ...
                dataFileName, err.message);
    end
end   % end of ff loop for batch operation, which started L110
% END OF BATCH LOOP
% 
% STEP 7 Write headers and matrix data into spread sheet.
CCC1 = {'xlsx = ',stackNameFF_noTif{1},'     ','     ','    ','      '};            % row1
CCC2 = {' ',' ',' ',' ','     ','      '};
CCC3 = {};
for ff = 1:numStacks
    CCC3 = [CCC3,stackNameFF{ff},' ',' ',' ','     ','     '];     
end
CCC4 = {};
for ff = 1:numStacks
    CCC4 = [CCC4,'   rr',' D',' dD',' dD/D ','     ','     '];           % row4
end
CCC5 = {};
for ff = 1:numStacks
    CCC5 = [CCC5,'    ',' nm^2/s',' nm^2/s',' ','     ','     ']; 
end
writecell(CCC1,fullFileNameOutputExcel,'Sheet',1,'Range','A1');
writecell(CCC2,fullFileNameOutputExcel,'Sheet',1,'Range','A2');
writecell(CCC3,fullFileNameOutputExcel,'Sheet',1,'Range','A3');
writecell(CCC4,fullFileNameOutputExcel,'Sheet',1,'Range','A4');
writecell(CCC5,fullFileNameOutputExcel,'Sheet',1,'Range','A5');
writematrix(excel_matrix, fullFileNameOutputExcel,'Sheet',1,'Range','A6');
toc;
end
% End of Main_Batch_7x7_V3.8