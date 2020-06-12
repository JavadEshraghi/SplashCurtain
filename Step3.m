clear all; clc; close all;
cd('F:\Processing\Sphere_experiment\Data_set_processing\Experimental_trajectory_info\375')
%% Directory
Sphere_Speed = '20'; % Impact velocity
Diameter = '0375'; % Sphere Diameter
Material = '1'; % the height it was dropped at
Trial = 2; % the number of the trial to be examined

folderName = ['S' num2str(Sphere_Speed) '_D' num2str(Diameter) '_M' num2str(Material)  '_0' num2str(Trial)];
load([folderName '_ExpTrajectory' '.mat'])
%% Initialization
% Splash Parameters
R = 0.00238125;
a = 0.14 * R;
Cd = 0.001;
deltaP = 150;
left = 1;
right = 0;
save_model_trajectory_right = 0;
save_model_trajectory_left = 0;

save_modified_trajectory = 0;
% Initial conditions
u_0 = 1.4;
v_0 = 1.7;

% Show comparing plots of two models (with and without 2nd surface tension)
show_comp_plot = 0;

% fluid properties and Constants
g = 9.81;
sigma = 0.0728;
rho = 998.21;
nu_air = 15.11e-6;
rhoair = 1.205;

%% Considering Pressure with constant Cd and delta P
f1 = @(t,x) [x(3);x(4);-(g*(x(3)*x(4))/(2*(x(3)^2+x(4)^2))+2*sigma*x(3)/(rho*pi*a^2*sqrt(x(3)^2+x(4)^2))+Cd*x(3)*sqrt(x(3)^2+x(4)^2)/(pi*a)+(deltaP*x(4))/(rho*pi*a*sqrt(x(3)^2+x(4)^2)));-(g*(x(3)^2+2*x(4)^2)/(2*(x(3)^2+x(4)^2))+2*sigma*x(4)/(rho*pi*a^2*sqrt(x(3)^2+x(4)^2))+Cd*x(4)*sqrt(x(3)^2+x(4)^2)/(pi*a)-(deltaP*x(3))/(rho*pi*a*sqrt(x(3)^2+x(4)^2)))];
[t1,xa1] = ode45(f1,[0:10e-5:0.0086],[R 0 u_0 v_0]);
r1 = xa1(:,1)/R;
z1 = xa1(:,2)/R;

f2 = @(t,x) [x(3);x(4);-(g*(x(3)*x(4))/(2*(x(3)^2+x(4)^2))+sigma*x(4)*(sqrt(x(3)^2+x(4)^2))*(2*x(3)^2+x(4)^2)/(2*a*x(1)*rho*pi*(x(3)^2+x(4)^2)^2)+2*sigma*x(3)/(rho*pi*a^2*sqrt(x(3)^2+x(4)^2))+Cd*sqrt(x(3)^2+2*x(4)^2)*x(3)/(pi*a)+(deltaP*x(4))/(rho*pi*a*sqrt(x(3)^2+x(4)^2)));-(g*(x(3)^2+2*x(4)^2)/(2*(x(3)^2+x(4)^2))+sigma*x(4)*(sqrt(x(3)^2+x(4)^2))*(x(3)*x(4))/(2*a*x(1)*rho*pi*(x(3)^2+x(4)^2)^2)+2*sigma*x(4)/(rho*pi*a^2*sqrt(x(3)^2+x(4)^2))+Cd*sqrt(x(3)^2+2*x(4)^2)*x(4)/(pi*a)-(deltaP*x(3))/(rho*pi*a*sqrt(x(3)^2+x(4)^2)))];
[t2,xa2] = ode45(f2,[0:10e-5:0.0086],[R 0 u_0 v_0]);
r2 = xa2(:,1)/R;
z2 = xa2(:,2)/R;

%% Plotting
if show_comp_plot == 1 
    figure('units','normalized','outerposition',[0 0 0.75 0.75])
    set(gca,'fontsize',14,'FontName','Garamond','FontWeight','bold','Color','w');
    set(gcf,'color','white');
    hold on
    plot(r1, z1, 'co-')
    hold on
    plot(r2, z2, 'ro-')
    xlabel('r/R','fontsize',16,'FontName','Garamond','FontWeight','bold')
    ylabel('z/R','fontsize',16,'FontName','Garamond','FontWeight','bold')
    legend('Neglecting 2nd surface tension','Considering 2nd surface tension','location','southeast')
    legend boxoff
    hold off
end

if left == 1
    figure('units','normalized','outerposition',[0 0 0.75 0.75]) 
    set(gca,'fontsize',14,'FontName','Garamond','FontWeight','bold','Color','w');
    set(gcf,'color','white');
    hold on
    plot(r2, z2, 'b-')
    hold on
    plot(rr_left_dimensionless,zz_left_dimensionless, 'rX')
    xlabel('r/R','fontsize',16,'FontName','Garamond','FontWeight','bold')
    ylabel('z/R','fontsize',16,'FontName','Garamond','FontWeight','bold')
    legend('Model Prediction','Experiment','location','southeast')
    legend boxoff
    hold off
end
if right == 1
    figure('units','normalized','outerposition',[0 0 0.75 0.75]) 
    set(gca,'fontsize',14,'FontName','Garamond','FontWeight','bold','Color','w');
    set(gcf,'color','white');
    hold on
    plot(r2, z2, 'b-')
    hold on
    plot(rr_right_dimensionless,zz_right_dimensionless, 'rX')
    xlabel('r/R','fontsize',16,'FontName','Garamond','FontWeight','bold')
    ylabel('z/R','fontsize',16,'FontName','Garamond','FontWeight','bold')
    legend('Model Prediction','Experiment','location','southeast')
    legend boxoff
    hold off
end

% Saving

if save_model_trajectory_right == 1
    % Save Info
    cd('F:\Processing\Sphere_experiment\Data_set_processing\Model_trajectory_info')
    save([folderName '_RightModelTrajectoryinfo' '.mat'])
end

if save_model_trajectory_left == 1
    % Save Info
    cd('F:\Processing\Sphere_experiment\Data_set_processing\Model_trajectory_info')
    save([folderName '_LeftModelTrajectoryinfo' '.mat'])
end


