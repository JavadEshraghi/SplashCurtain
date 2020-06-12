%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% either bofore or after of this step we need to determine the sphere radius
% and x position of the center of the sphere to be able to dimensionalize
% and find the mapping function to compare the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
clc; 
clear;

%% Directory
Sphere_Speed = '35'; % Impact velocity
Diameter = '0750'; % Sphere Diameter
Material = '2'; % the height it was dropped at
Trial = 3; % the number of the trial to be examined
Sphere_Diameter_Splash_Im = 180; % Sphere diameter in Pixel

ssFrame = 22; % the frame where the sphere completely appeared
ccFrame = 88; % the frame that the sphere contacted the surface
closureFrame = 265; % the frame that the splash is closed

imageName_submerged = sprintf('S%s_D%s_M%s_%02d%05d.tif',Sphere_Speed,Diameter,Material,Trial,ssFrame);
imageName_contact = sprintf('S%s_D%s_M%s_%02d%05d.tif',Sphere_Speed,Diameter,Material,Trial,ccFrame);
folderName = ['S' num2str(Sphere_Speed) '_D' num2str(Diameter) '_M' num2str(Material)  '_0' num2str(Trial)];
cd(num2str(folderName))

figure(1)
imshow(imageName_submerged)
[x_right_side, y_right_side] = ginput(1);
Sphere_right_side = x_right_side;

figure(2)
imshow(imageName_contact)
[x_level, y_level] = ginput(1);



Sphere_Radius = Sphere_Diameter_Splash_Im/2;
Sphere_Center_X = Sphere_right_side - Sphere_Radius;
Im_Size = size(imageName_submerged);
Image_height = Im_Size(2);
Free_Surface_Y = y_level;
set_height = 1;

% cd('F:\..\Sphere_experiment\Data_set_images\Splash_camera')


%% 
calib = Sphere_Diameter_Splash_Im/(str2double(Diameter)*2.54/1000); %pix/cm
rho_air = 1.205;

% %% Initializations
% Cavity_volume = zeros;

%% Select the surface height
% figure(3)
% Im1 = imread(imageName_contact);
% imshow(Im1)
% [x1, y1] = ginput(1);

% %% Main Time Loop, Frame by Frame increments
% for ii = ssFrame:ccFrame
% 
%     Im_Name = sprintf('S%s_D%s_M%s_%02d%06d.tif',Sphere_Speed,Diameter,Material,Trial,ii);
%     Im = imread(Im_Name); % reads the image
%     Im_Size_orig = size(Im);
%     
%     %% Image Processing
%     BWcrop = imcrop(Im,[1 1 Im_Size_orig(2) y1]); % crops the image
%     BW = im2bw(BWcrop,.7); % converts image to black and white
%     % only pixels that are connected to another # pixels are kept
%     BWred = bwareaopen(~BW, 100);
%     BWfill=imfill(BWred,'holes'); % filles in holes of the image   
%     B = bwboundaries(BWfill); % finds the boundaries of the image
%     sizeB = size(B); % size of the boundary matrix
%     
%     %%%% Sometimes, B will have empty cell array (sizeB = [0 1])
%     if sizeB ~= 0
%         bound = B{1};
%     else
%         bound = [0 0];
%     end
% 
%     % Stores the boundary in variables
%     boundY = bound(:,1)';
%     boundX = bound(:,2)';
%     
%     % Find velocity of the spheres
%     S_y_max (ii) = max(boundY);
% end
% 
% Sphere_Impact_Vel = diff (S_y_max); % px per frame
% Sphere_Impact_Velocity = Sphere_Impact_Vel(end);


%% Find Trajectory and initial velocity

% rr_right_dimensionless = (rim_x_right - Sphere_Center_X) / Sphere_Radius;
rr_left_dimensionless = -(rim_x_left - Sphere_Center_X) / Sphere_Radius;
% zz_right_dimensionless = ((Free_Surface_Y) - rim_y_right) / Sphere_Radius ;
zz_left_dimensionless = ((Free_Surface_Y) - rim_y_left) / Sphere_Radius ;


cd ('F:\Processing\Sphere_experiment\Data_set_processing\Experimental_trajectory_info')
save([folderName '_ExpTrajectory' '.mat'])

