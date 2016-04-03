% given the length of the station array L_s and the average Moho
% depth H, figure out the appropriate geometry for a plane wave
% with incident angle th (in degrees) that will propagate (as
% initial condition to the free surface and cover the entire array

close all; clear all;
L_s=400; %km
H=40; %km
L_1=1./4*L_s; % km
angle_inc=23;
th=angle_inc*pi/180;

% plane waves start below the Moho
L_1=max(L_1,H*tan(th));

x0=0.; z0=0.;
x1=x0+L_1;
x2=x1+L_s;
xm=x2+x1;

% 
fac=1.2;
xp=x0-H*fac/tan(th);
x4=x1-(x1-xp)*(sin(th).^2);
z4=(x1-xp)*sin(th)*cos(th);
x5=x2-(x2-xp)*(sin(th).^2);
z5=(x2-xp)*sin(th)*cos(th);
zm=z5*1.1;

% draw the setup geometry
set(0,'Defaultlinelinewidth',2);
rectangle('position',[x0,z0,xm,zm],'EdgeColor','g'); hold on;
% extending -x axis
plot([x0,xp],[zm,zm],'g:');
% rays going to x1 and x2 receivers
plot([x1,x4],[zm,zm-z4]);
plot([x2,x5],[zm,zm-z5]);
% wavefront connecting (x4,z4), (x5,z5)
plot([xp,x4,x5],[zm,zm-z4,zm-z5],'r:');
% plot Moho
plot([x0,xm],[zm-H,zm-H],'m')
axis equal;

xp,xm,zm