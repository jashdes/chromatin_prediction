function [Xc, Yc, theta_mean, theta_mean_std] = center_of_rotation_func(x_cor, y_cor)  %x_cor and y_cor are the 9x200 matrix with all drift-corrected x and y values
% center_of_rotation_func.m
% gholz March 2017
% based on Matlab Central 5332l9 Yong-San Yoon, KAIST for the 3D case
% modified for the 2D case by gholzwarth

%  Assume m points in each frame, n frames, 2 dimensions x,y
%  Motivation: find an approach which sets up a linear equation Ax = b, with x
%  the unknown coordinates Xc and Yc for the center of rotation.
%  [Xc, Yc]: center of rotation
%  xji, yji: Cartesian coordinates of the j.th point in the i.th frame.
%  Let rji^2 = squared distance from point j in frame i to center of
%  rotation. If Xc and Yc are correct, then
%       r11^2 = r12^2;
%       r21^2 = r22^2;
% writing out the resulting details in matrix form, Yoon gets

%  / 2*x11 2*y11 1 0 ... \          / x11^2 + y11^2 \
%  | 2*x12 2*y12 1 0 ... | / Xc \   | x12^2 + y12^2 |
%  | 2*x13 2*y13 1 0 ... | | Yc |   | x13^2 + y13^2 |
%  | 2*x1n 2*y1n 1 0 ... | | P1 |   | x1n^2 + y1n^2|
%  |                     | | P2 |   | ...                   |
%  | 2*x21 2*y21 0 1 ... | | .. | = | x21^2 + y21^2|
%  | 2*x22 2*y22 0 1 ... | | .. |   | x22^2 + y22^2|
%  | 2*x23 2*y23 0 1 ... | | .. |   | x23^2 + y23^2|
%  | ...                       | | .. |   | ...                   |
%  | 2*xm1 2*ym1 ... 0 1 | | .. |   | xm1^2 + ym1^2|
%  | ...                       | \ Pm /   | ...                   |
%  \ 2*xmn 2*ymn ... 0 1 /          \ xmn^2 + ymn^2/

%  where 
%  P<j>: Residuals of radii: P<j> = R<j>^2 - (Xc^2 + Yc^2)

% image 1, 9 spots
% define the x- and y-data for the original line we would like to rotate
%     x = [-.5 -.5 -.5 0.0 0.0 0.0 0.5 0.5 0.5];
%     y = [0.5 0.0 -.5 0.5 0.0 -.5 0.5 0.0 -.5];  
% create a matrix of these points, which will be useful in future calculations
%     v = [x;y];
% choose a point which will be the center of rotation
% x_center = 0.765;
% y_center = 0.235;

% create a matrix which will be used later in calculations
%     center = repmat([x_center; y_center], 1, length(x));
% set rotation angle  (deg)    
%     theta = 10; % degrees
% construct rotation matrix
%     R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)]
% Apply rotation about center and restore to original coordinate system 
%     vo = R*(v - center) + center
% pick out the vectors of rotated x- and y-data
% x_rotated = vo(1,:);
% y_rotated = vo(2,:);
% make a plot
% plot(x, y, 'k-', x_rotated , y_rotated, 'r-', x_center, y_center, 'bo');
% axis equal
% define p1:p9 in format for matrix inversion
p1 = [x_cor(1,1), y_cor(1,1); x_cor(1,50), y_cor(1,50)];
p2 = [x_cor(2,1), y_cor(2,1); x_cor(2,50), y_cor(2,50)];
p3 = [x_cor(3,1), y_cor(3,1); x_cor(3,50), y_cor(3,50)];
p4 = [x_cor(4,1), y_cor(4,1); x_cor(4,50), y_cor(4,50)];
p5 = [x_cor(5,1), y_cor(5,1); x_cor(5,50), y_cor(5,50)];
p6 = [x_cor(6,1), y_cor(6,1); x_cor(6,50), y_cor(6,50)];
p7 = [x_cor(7,1), y_cor(7,1); x_cor(7,50), y_cor(7,50)];
p8 = [x_cor(8,1), y_cor(8,1); x_cor(8,50), y_cor(8,50)];
p9 = [x_cor(9,1), y_cor(9,1); x_cor(9,50), y_cor(9,50)];
%pause ()

%%Plot points
% figure (1)
% % plot 9 points in image 1
% plot(p1(1,1),p1(1,2),'or');hold on;
% plot(p2(1,1),p2(1,2),'og');hold on;
% plot(p3(1,1),p3(1,2),'ob');hold on;
% plot(p4(1,1),p4(1,2),'or');hold on;
% plot(p5(1,1),p5(1,2),'or');hold on;
% plot(p6(1,1),p6(1,2),'or');hold on;
% plot(p7(1,1),p7(1,2),'or');hold on;
% plot(p8(1,1),p8(1,2),'or');hold on;
% plot(p9(1,1),p9(1,2),'or');hold on;
% % plot 9 points in image 2
% plot(p1(2,1),p1(2,2),'sqr');hold on;
% plot(p2(2,1),p2(2,2),'sqb');hold on;
% plot(p3(2,1),p3(2,2),'sqg');hold on;
% plot(p4(2,1),p4(2,2),'sqb');hold on;
% plot(p5(2,1),p5(2,2),'sqb');hold on;
% plot(p6(2,1),p6(2,2),'sqb');hold on;
% plot(p7(2,1),p7(2,2),'sqb');hold on;
% plot(p8(2,1),p8(2,2),'sqb');hold on;
% plot(p9(2,1),p9(2,2),'sqb');hold on;
% title ('Fig1 9 pts, 2 frames ')

% assume 2D case, 2 frames, 2 points.  Check the Yoon method.
        % %point 1
        %     p1 = [0, 1 ; 1 , 0  ] % notation: [p1(x,image1) p1(y(image1) ; p1(x,image2) p1(y,image2)]
        % % point 2
        %     p2 = [1, 0   ; 0,  -1]
%% Adding noise
% noise1 = 0.1*randn(2,2)
% noise2 = 0.1*randn(2,2) 
% noise3 = 0.1*randn(2,2) 
% noise4 = 0.1*randn(2,2) 
% noise5 = 0.1*randn(2,2) 
% noise6 = 0.1*randn(2,2) 
% noise7 = 0.1*randn(2,2) 
% noise8 = 0.1*randn(2,2) 
% noise9 = 0.1*randn(2,2) 
% 
% p1 = p1+noise1
% p2 = p2+noise2
% p3 = p3+noise3
% p4 = p4+noise4
% p5 = p5+noise5
% p6 = p6+noise6
% p7 = p7+noise7
% p8 = p8+noise8
% p9 = p1+noise9
% pause ()

%  p1 = [0   0.5 ; 0.5   0  ]
%  p2 = [0.5 0   ; 0    -0.5]
% point 1
%    p1 = [0, 1 ; 1 , 0  ]; % notation: [p1(x,image1) p1(y(image1) ; p1(x,image2) p1(y,image2)]
% point 2
%    p2 = [0, 2  ; 2,  0];


    %     p1 = rand(2, 2);  % [2 frames, 2 coordinates, 4 points ]
    %     p2 = rand(2, 2);
    % p3 = rand(2, 2);
    % p4 = rand(2, 2);
nPoint = 9;
nFrame = 2;
MPM = cat(1, p1, p2, p3, p4, p5, p6, p7, p8, p9);      % multiple points as matrix

% Define matrices A and b such that Ax   = b;
A   = cat(2, 2 .* MPM, zeros(nFrame * nPoint, nPoint));
for iN = 1:nPoint
  A(((iN - 1) * nFrame + 1):(iN * nFrame), 3 + iN) = 1;
end
b   = sum(MPM .* MPM, 2);

% solve for x formally
x  = A \ b;  % This shows the power of Matlab to solve this problem efficiently
C  = transpose(x(1:3));  % C(1) = Xcenter, C(2) = Ycenter.
Xc = C(1);
Yc = C(2);
        % fprintf('display Xc and Yc');
        % display(Xc);
        % display(Yc);

Res = A * x - b;

% find angle of rotation using dot product formula for vectors
    % v1 dot v2 = v1*v2*cos(theta)
% for spot 1, find vector v1 from Xc,Yc to p1 in frame 1
     v1 = [Xc, Yc] - [p1(1,1),p1(1,2)];  % frame 1
% for spot 1, find vector v2 from Xc,Yc to p1 in frame 2     
     v2 = [Xc, Yc] - [p1(2,1),p1(2,2)];  % frame 2
% compute angle of rotation from the dot product of v1 and v2 both normalized
     theta1 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
             % fprintf('theta1(degrees) = ');
             % display(theta1);
% repeat for spot 2:
     v1 = [Xc, Yc] - [p2(1,1),p2(1,2)];
     v2 = [Xc, Yc] - [p2(2,1),p2(2,2)];
     theta2 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
            % fprintf('theta2(degrees) = ');
            % display(theta2);
     
% repeat for spot 3:
     v1 = [Xc, Yc] - [p3(1,1),p3(1,2)];
     v2 = [Xc, Yc] - [p3(2,1),p3(2,2)];
     theta3 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
        %      fprintf('theta3(degrees) = ');
        %      display(theta3);
        %      
% repeat for spot 4:
     v1 = [Xc, Yc] - [p4(1,1),p4(1,2)];
     v2 = [Xc, Yc] - [p4(2,1),p4(2,2)];
     theta4 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
        % fprintf('theta4(degrees) = ');
        % display(theta4);
     
% repeat for spot 5:
     v1 = [Xc, Yc] - [p5(1,1),p5(1,2)];
     v2 = [Xc, Yc] - [p5(2,1),p5(2,2)];
     theta5 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
            % fprintf('theta5(degrees) = ');
            % display(theta5);
     
% repeat for spot 6:
     v1 = [Xc, Yc] - [p6(1,1),p6(1,2)];
     v2 = [Xc, Yc] - [p6(2,1),p6(2,2)];
     theta6 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
            % fprintf('theta6(degrees) = ');
            % display(theta6);
     
% repeat for spot 7:
     v1 = [Xc, Yc] - [p7(1,1),p7(1,2)];
     v2 = [Xc, Yc] - [p7(2,1),p7(2,2)];
     theta7 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
            % fprintf('theta7(degrees) = ');
            % display(theta7);
     
% repeat for spot 8:
     v1 = [Xc, Yc] - [p8(1,1),p8(1,2)];
     v2 = [Xc, Yc] - [p8(2,1),p8(2,2)];
     theta8 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
           % fprintf('theta8(degrees) = ');
           % display(theta8);
     
% repeat for spot 9:
     v1 = [Xc, Yc] - [p9(1,1),p9(1,2)];
     v2 = [Xc, Yc] - [p9(2,1),p9(2,2)];
     theta9 = acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
            % fprintf('theta9(degrees) = ');
            % display(theta9);
% evaluate theta_mean and theta_mean_std
     theta_array = [theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8,theta9];
     theta_mean  = mean(theta_array);
     theta_mean_std = std(theta_array);
   

% The more points, the better the results. But noise and translation will be a severe problem: If you mark a bunch of points on a rigid body, rotate it around several axes (or around a point in 2D). If there is a translation in addition, the found "center of rotation" depends on the position of the marked points...