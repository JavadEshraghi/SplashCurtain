close all; 
clc; 
clear;

%% Directory
Sphere_Speed = '20'; % Impact velocity
Diameter = '0750'; % Sphere Diameter
Material = '4'; % density ratio
Trial = 1; % trial #
Sphere_Diameter_Cavity_Im = 85; % Sphere diameter in Pixel
fps = 5000; % the rate that the camera was filming at
Thresh = 50000;

cFrame = 10; % the frame where contact with the surface occurs
sFrame = 15; % the frame where the sphere completely submerged
pFrame = 334; % the frame that pinch off occurs at
derivative_order = 4; % The order of numerical derivative to find air flow rate

imageName_contact = sprintf('S%s_D%s_M%s_%02d%05d.tif',Sphere_Speed,Diameter,Material,Trial,cFrame);
imageName_pinchoff = sprintf('S%s_D%s_M%s_%02d%05d.tif',Sphere_Speed,Diameter,Material,Trial,pFrame);

folderName = ['S' num2str(Sphere_Speed) '_D' num2str(Diameter) '_M' num2str(Material)  '_0' num2str(Trial)];

% cd('F:\..\Sphere_experiment\Data_set_images\Cavity_camera')
% cd(num2str(folderName))


%% Body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User selects the following points on the submerged (one point) and
% pinch-off (three points):
% i.	free surface
% ii.	pinch-off location
% iii.	bottom of sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
Im1 = imread(imageName_contact);
imshow(Im1)
[x1, y1] = ginput(1);

figure(2)
Im2 = imread(imageName_pinchoff);
imshow(Im2)
[x2, y2] = ginput(2);

freesurface_Y = y1;
sphere_center_Y = (y2(2)-Sphere_Diameter_Cavity_Im/2);
pinchoff_Y = y2(1);


%% Initializations
calib = Sphere_Diameter_Cavity_Im/(0.75*2.54); %pix/cm
R = Sphere_Diameter_Cavity_Im/2; % radius, pi
rho_air = 1.205;
Cavity_volume = zeros;

%% Main Time Loop, Frame by Frame increments
for ii = sFrame:pFrame
    Im_Name = sprintf('S%s_D%s_M%s_%02d%05d.tif',Sphere_Speed,Diameter,Material,Trial,ii);
    Im = imread(Im_Name); % reads the image
    Im_Size_orig = size(Im);
    %% Image Processing
    BWcrop = imcrop(Im,[1 y1 Im_Size_orig(2) Im_Size_orig(1)]); % crops the image 
    %%Volume calculation
    left = zeros;
    right = zeros;
    Cavity_radius = zeros;
    Im_Size = size(BWcrop);
    
    for jj = 1:Im_Size(1)
        for kk = 1:Im_Size(2)
            if BWcrop(jj,kk)>Thresh
                left (jj) = kk;
                break
            end
        end
        for kk = 1:Im_Size(2)            
            if BWcrop(jj,Im_Size(2)-kk+1)>Thresh
                right (jj) = Im_Size(2)-kk;
                break
            end
        end
    end
    % Find velocity of the spheres
    y_max (ii) = max(length(left),length(right)); 
    
    imshow(BWcrop)
    hold on
    plot(left,'LineWidth',2)
    plot(right,'LineWidth',2)
    hold off
    pause(0.005)
    Cavity_radius = (right-left)/2;
    Cavity_volume (ii) = sum(pi.*Cavity_radius.^2)-(4/3*pi*R^3);   
end
Vol_Smoothed = smooth(0:1:(pFrame-1),Cavity_volume(1:pFrame),0.35,'rloess'); % Smoothing
Vol_Smoothed_P = Vol_Smoothed;
% if Vol_Smoothed(1)<0
   Vol_Smoothed = Vol_Smoothed - min(Vol_Smoothed);
% end
Volume_Dimensionless = Vol_Smoothed/(4/3*pi*R^3); 

% 2nd order derivative
if (derivative_order==2)
    for jj = 2:(pFrame-1)
        Vol_Flow_Rate (jj) = 1/2 * (Vol_Smoothed(jj+1) - Vol_Smoothed(jj-1)); % cubic px per frame
    end
end

% 4th order derivative
if (derivative_order==4)
    for jj = 2:(pFrame-3)
        Vol_Flow_Rate (jj-1) = 1/12 * (-3 * Vol_Smoothed(jj-1) - 10 * Vol_Smoothed(jj) + 18 * Vol_Smoothed(jj+1) - 6 * Vol_Smoothed(jj+2) + Vol_Smoothed(jj+3));
    end
end

Vol_Flow_Rate_size=size(Vol_Flow_Rate);
Vol_flow_rate_Smoothed = smooth(1:1:(Vol_Flow_Rate_size(2)),Vol_Flow_Rate(1:Vol_Flow_Rate_size(2)),0.35,'rloess'); % Smoothing
Volume_flow_rate = Vol_flow_rate_Smoothed; % cubic px per frame
Air_velocity = Volume_flow_rate/(pi*R^2); % px per frame
Air_velocity_Dimensionless = Air_velocity / ((str2double(Sphere_Speed)*calib/fps*10));
Pressure_diff = rho_air * 0.5 * Air_velocity.^2; % kg/(px.frame squared) 
Pressure_diff_Dimensionless = Pressure_diff / (rho_air * 0.5 * (str2double(Sphere_Speed)*calib/fps*10).^2); %Dimensionless
Sphere_Vel = diff (y_max); % px per frame
Sphere_Pinchoff_Vel = mean(Sphere_Vel(end-5:end));

figure,plot(Volume_Dimensionless)
figure,plot(Cavity_volume)
hold on
plot(Vol_Smoothed)
hold off
figure,plot(Air_velocity_Dimensionless)
figure,plot(Pressure_diff_Dimensionless)
% Save Functions
cd('F:\Processing\Sphere_experiment\Data_set_processing\Cavity_volume_info')
% save([folderName '_CavityVolumeinfo' '.mat'],'Sphere_Pinchoff_Vel','Sphere_Vel','Pressure_diff_Dimensionless','Cavity_volume','Volume_flow_rate','Air_velocity')
save([folderName '_Cavityinfo' '.mat'])
