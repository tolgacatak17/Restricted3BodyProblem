%=================================%
% ME 303 PROJECT ASSIGNMENT 2
% Question 1
% 2018405144-Tolga Çatak

% Restricted 3-Body Problem

clear all
close all
clc

format long

%%
%======Define the given parameters======%

% Define the mass of the Moon and the Earth
m1 = 0.012277471;
m2 = 1-m1;

% Define the time interval and the step size
T = 17.06521656015796;
tmin = 0;
tmax = T;

h_e = T/24000;
h_rk4 = T/6000;

t_e = tmin:h_e:tmax;
t_rk4 = tmin:h_rk4:tmax;

% Define the number of iteration
M_rk4 = (tmax-tmin)/h_rk4;
M_e = (tmax-tmin)/h_e;

%%
%======Define the equations======%

% Define the equations as a system of first order differential equations

%z1=y1
%z2=y2
%z3=y1'
%z4=y2'

%z1'=y1'
%z2'=y2'
%z3'=y1''
%z4'=y2''

% Define the slope functions
f1 = @(t,z1,z2,z3,z4) z3;
f2 = @(t,z1,z2,z3,z4) z4;
f3 = @(t,z1,z2,z3,z4) z1 + 2*z4 - m2*(z1+m1)/((z1+m1)^2+z2^2)^(3/2) - m1*(z1-m2)/((z1-m2)^2+z2^2)^(3/2);
f4 = @(t,z1,z2,z3,z4) z2 - 2*z3 - m2*z2/((z1+m1)^2+z2^2)^(3/2) - m1*z2/((z1-m2)^2+z2^2)^(3/2);

% Preallocate space
z1 = zeros(1,M_e+1);
z2 = zeros(1,M_e+1);
z3 = zeros(1,M_e+1);
z4 = zeros(1,M_e+1);

% Define the initial values
z1_0 = 0.994;
z2_0 = 0;
z3_0 = 0;
z4_0 = -2.0015851063790825;

% Assign the initial values to the arrays
z1(1)=z1_0;
z2(1)=z2_0;
z3(1)=z3_0;
z4(1)=z4_0;

%%
%=======Euler Method======%

% Apply the Euler Method
for k=1:M_e+1
    z1(k+1)=z1(k)+h_e*f1(t_e(k),z1(k),z2(k),z3(k),z4(k));
    z2(k+1)=z2(k)+h_e*f2(t_e(k),z1(k),z2(k),z3(k),z4(k));
    z3(k+1)=z3(k)+h_e*f3(t_e(k),z1(k),z2(k),z3(k),z4(k));
    z4(k+1)=z4(k)+h_e*f4(t_e(k),z1(k),z2(k),z3(k),z4(k));
end

% Plot the results of the Euler Method
subplot(1,3,1)
hold on
plot(z1,z2,'b',LineWidth=2)
% Define the boundaries of the coordinate system
xmin=-2; xmax=2; ymin=-1.5; ymax=2.5;
axis([xmin xmax ymin ymax])
grid on
plot(0,0, '.k', 'linewidth', 4)
title('Euler Solution')
xlabel('y1')
ylabel('y2')
legend('Trajectory', 'Origin');

% Store the values for the animation.
z1e=z1;
z2e=z2;

%%
%======Runge-Kutta 4 Method======%

% Reset the values of the arrays to use the same variables
z1 = zeros(1,M_rk4+1);
z2 = zeros(1,M_rk4+1);
z3 = zeros(1,M_rk4+1);
z4 = zeros(1,M_rk4+1);

% Assign the initial values to the arrays
z1(1)=z1_0;
z2(1)=z2_0;
z3(1)=z3_0;
z4(1)=z4_0;

% Apply the Runge-Kutta of Order N=4 Method
for k=1:M_rk4+1
    f1_1=f1(t_rk4(k),z1(k),z2(k),z3(k),z4(k));
    f2_1=f2(t_rk4(k),z1(k),z2(k),z3(k),z4(k));
    f3_1=f3(t_rk4(k),z1(k),z2(k),z3(k),z4(k));
    f4_1=f4(t_rk4(k),z1(k),z2(k),z3(k),z4(k));
    
    f1_2=f1(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_1,z2(k)+h_rk4/2*f2_1,z3(k)+h_rk4/2*f3_1,z4(k)+h_rk4/2*f4_1);
    f2_2=f2(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_1,z2(k)+h_rk4/2*f2_1,z3(k)+h_rk4/2*f3_1,z4(k)+h_rk4/2*f4_1);
    f3_2=f3(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_1,z2(k)+h_rk4/2*f2_1,z3(k)+h_rk4/2*f3_1,z4(k)+h_rk4/2*f4_1);
    f4_2=f4(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_1,z2(k)+h_rk4/2*f2_1,z3(k)+h_rk4/2*f3_1,z4(k)+h_rk4/2*f4_1);
    
    f1_3=f1(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_2,z2(k)+h_rk4/2*f2_2,z3(k)+h_rk4/2*f3_2,z4(k)+h_rk4/2*f4_2);
    f2_3=f2(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_2,z2(k)+h_rk4/2*f2_2,z3(k)+h_rk4/2*f3_2,z4(k)+h_rk4/2*f4_2);
    f3_3=f3(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_2,z2(k)+h_rk4/2*f2_2,z3(k)+h_rk4/2*f3_2,z4(k)+h_rk4/2*f4_2);
    f4_3=f4(t_rk4(k)+h_rk4/2,z1(k)+h_rk4/2*f1_2,z2(k)+h_rk4/2*f2_2,z3(k)+h_rk4/2*f3_2,z4(k)+h_rk4/2*f4_2);

    f1_4=f1(t_rk4(k)+h_rk4,z1(k)+h_rk4*f1_3,z2(k)+h_rk4*f2_3,z3(k)+h_rk4*f3_3,z4(k)+h_rk4*f4_3);
    f2_4=f2(t_rk4(k)+h_rk4,z1(k)+h_rk4*f1_3,z2(k)+h_rk4*f2_3,z3(k)+h_rk4*f3_3,z4(k)+h_rk4*f4_3);
    f3_4=f3(t_rk4(k)+h_rk4,z1(k)+h_rk4*f1_3,z2(k)+h_rk4*f2_3,z3(k)+h_rk4*f3_3,z4(k)+h_rk4*f4_3);
    f4_4=f4(t_rk4(k)+h_rk4,z1(k)+h_rk4*f1_3,z2(k)+h_rk4*f2_3,z3(k)+h_rk4*f3_3,z4(k)+h_rk4*f4_3);

    % Calculate the approximated function values for each iteration
    z1(k+1)=z1(k)+h_rk4/6*(f1_1+2*f1_2+2*f1_3+f1_4);
    z2(k+1)=z2(k)+h_rk4/6*(f2_1+2*f2_2+2*f2_3+f2_4);
    z3(k+1)=z3(k)+h_rk4/6*(f3_1+2*f3_2+2*f3_3+f3_4);
    z4(k+1)=z4(k)+h_rk4/6*(f4_1+2*f4_2+2*f4_3+f4_4);
end

% Plot the results of the RK4 Method
subplot(1,3,2)
hold on
plot(z1,z2,'c',LineWidth=2)
axis([xmin xmax ymin ymax])
grid on
plot(0,0, '.k', 'linewidth', 4)
title('RK4 Solution')
xlabel('y1')
ylabel('y2')
legend('Trajectory', 'Origin');

%%
%======ODE45 Method======%

% Define a function array for ode45 function
f= @(t,z) [z(3);z(4);z(1)+2*z(4)-m2*(z(1)+m1)/((z(1)+m1)^2+z(2)^2)^(3/2)-m1*(z(1)-m2)/((z(1)-m2)^2+z(2)^2)^(3/2);z(2)-2*z(3)-m2*z(2)/((z(1)+m1)^2+z(2)^2)^(3/2)-m1*z(2)/((z(1)-m2)^2+z(2)^2)^(3/2)];
% Define an initial value array
z0 = [z1_0;z2_0;z3_0;z4_0];
% Use the ode45 function
[t_ode45,z_ode45]=ode45(f,[tmin tmax],z0);

% Plot the results of the ode45 Method
subplot(1,3,3)
hold on
plot(z_ode45(:,1),z_ode45(:,2),'m',LineWidth=2)
axis([xmin xmax ymin ymax])
grid on
plot(0,0, '.k', 'linewidth', 4)
title('ode45 Solution')
xlabel('y1')
ylabel('y2')
legend('Trajectory', 'Origin');

% The plots are shown on the same figure with the same coordinate system
% Thus, the trajectories are to scale and can be compared

%%
%======Euler Method Animation======%

% Extract every nth element of the arrays to arrange the frames.
% This part is necessary since we are using RK4 method to solve the problem
% Normally, frame rate can be adjusted but it works by changing the step size of the solution 
% We apply this method to keep the step size of the RK4 Method constant while being able to change the frame rate of the animation
% As an increases, the framerate decreases, i.e. the resolution decreases

n=100;
z1e = z1e(1 : n : end);  
z2e = z2e(1 : n : end);  

% Use cool function to take the colormap
c = cool(6); 

figure
% Adjust the pixels of the animation
set(gcf,'Position',[50 50 640 640])     % Social

% Adjust the coordinate system of the animation
hold on ; grid on ; box on ; axis equal
set(gca,'XLim',[xmin xmax])
set(gca,'YLim',[ymin ymax])
set(gca,'XTick',[],'YTick',[])

%Choose a color for the background
set(gca,'Color','w')

% Create and open video writer object
v = VideoWriter('3body_euler.mp4','MPEG-4');
v.Quality   = 10;

% Start the animation
open(v);

grid on
title('Euler Solution Animation')
xlabel('y1')
ylabel('y2')

for i = 1:length(z1e)
    cla

    % Animate each of the components
    % Choose any desired color for the trajectory
    
    % Origin
    p = plot(0,0,'wo');
    set(p,'MarkerFaceColor','r','MarkerSize',5)
    
    % Trajectory Euler
    plot(z1e(1:i),z2e(1:i),'color',c(3,:),'LineWidth',3)
    
    % Position Euler
    p = plot(z1e(i),z2e(i),'r');
    set(p,'Marker','o','MarkerFaceColor',c(3,:),'Color','k','MarkerSize',5)
   
    frame = getframe(gcf);
    writeVideo(v,frame);
end

% Close the animation
close(v);

%%
%======RK4 Method Animation======%

% This part is used for the same reason, i.e. for frame rate.
n=100;
z1 = z1(1 : n : end);  
z2 = z2(1 : n : end);  

% Use cool function to take the colormap
c = cool(6); 

figure
% Adjust the pixels of the animation
set(gcf,'Position',[50 50 640 640])     % Social

% Adjust the coordinate system of the animation
hold on ; grid on ; box on ; axis equal
set(gca,'XLim',[xmin xmax])
set(gca,'YLim',[ymin ymax])
set(gca,'XTick',[],'YTick',[])

%Choose a color for the background
set(gca,'Color','w')

% Create and open video writer object
v = VideoWriter('3body_rk4.mp4','MPEG-4');
v.Quality  = 10;

% Start the animation
open(v);
    
title('RK4 Solution Animation')
xlabel('y1')
ylabel('y2')

for i = 1:length(z1)
    cla

    % Animate each of the components
    % Choose any desired color for the trajectory
    
    % Origin
    p = plot(0,0,'wo');
    set(p,'MarkerFaceColor','r','MarkerSize',5)
    
    % Trajectory RK4
    plot(z1(1:i),z2(1:i),'color',c(1,:),'LineWidth',3)
    
    % Position RK4
    p = plot(z1(i),z2(i),'r');
    set(p,'Marker','o','MarkerFaceColor',c(1,:),'Color','k','MarkerSize',5)
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end

% Close the animation
close(v);

%%
%======ode45 Method Animation======%

z1_ode45 = z_ode45(:,1);
z2_ode45 = z_ode45(:,2);

% This part is used for the same reason, i.e. for frame rate.
n=2;
z1_ode45 = z1_ode45(1 : n : end);  
z2_ode45 = z2_ode45(1 : n : end);  

% Use cool function to take the colormap
c = cool(6); 

figure
% Adjust the pixels of the animation
set(gcf,'Position',[50 50 640 640])     % Social

% Adjust the coordinate system of the animation
hold on ; grid on ; box on ; axis equal
set(gca,'XLim',[xmin xmax])
set(gca,'YLim',[ymin ymax])
set(gca,'XTick',[],'YTick',[])

%Choose a color for the background
set(gca,'Color','w')

% Create and open video writer object
v = VideoWriter('3body_ode45.mp4','MPEG-4');
v.Quality  = 10;

% Start the animation
open(v);

title('ode45 Solution Animation')
xlabel('y1')
ylabel('y2')

for i = 1:length(z1_ode45)
    cla

    % Animate each of the components
    % Choose any desired color for the trajectory
    
    % Origin
    p = plot(0,0,'wo');
    set(p,'MarkerFaceColor','r','MarkerSize',5)
    
    % Trajectory ode45
    plot(z1_ode45(1:i),z2_ode45(1:i),'color',c(6,:),'LineWidth',3)
    
    % Position ode45
    p = plot(z1_ode45(i),z2_ode45(i),'r');
    set(p,'Marker','o','MarkerFaceColor',c(6,:),'Color','k','MarkerSize',5)
    
    title('ode45 Solution Animation')
    xlabel('y1')
    ylabel('y2')

    frame = getframe(gcf);
    writeVideo(v,frame);
end

% Close the animation
close(v);