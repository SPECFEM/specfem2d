clear all; close all;
LENGTH=250000; L1=LENGTH/1000; H=40;
angle_inc=18;
th=angle_inc*pi/180;
dirr='OUTPUT_FILES_FORWARD/';
dir = 'OUTPUT_FILES_KERNEL/';
% receiver
file=[dirr 'receiver.dat'];
rec=load(file);
rec=rec/LENGTH; nrec=length(rec);

% kernel file
file=[dir 'proc000000_rhop_alpha_beta_kernel.dat'];
s=load(file);
xx=s(:,1)/LENGTH;zz=s(:,2)/LENGTH;
minx=min(xx); maxx=max(xx); minz=min(zz); maxz=max(zz);
nx=400; nz=400;
% define mesh
xa=linspace(minx,maxx,nx);za=linspace(minz,maxz,nz);
[X,Z] = meshgrid(xa,za);
rhop=s(:,3); alpha=s(:,4); beta=s(:,5);

subplot(3,1,1)
%scatter(xx,zz,30,rhop,'filled'); hold on;
FV=griddata(xx,zz,rhop,X,Z);
imagesc(xa,za,FV); set(gca,'YDir','normal'); hold on;
plot(rec(:,1),rec(:,2),'kx')
text(0.2,0.2,'K_{\rho_p} ')
quiver(0.2,1-H/L1, 0.1,0,'MaxHeadsize', 30)
quiver(0,0,0.4*sin(th),0.4*cos(th),'MaxHeadsize', 30)
title(['inc angle=', num2str(angle_inc), ';  nrec=', num2str(nrec)]);
jpeg_header='right-kernel-';
jpeg_file=[jpeg_header num2str(angle_inc) '-' num2str(nrec) '.jpg'];

subplot(3,1,2)
%scatter(xx,zz,30,beta,'filled'); hold on;
FV=griddata(xx,zz,beta,X,Z);
imagesc(xa,za,FV); set(gca,'YDir','normal'); hold on;
plot(rec(:,1),rec(:,2),'kx')
plot(rec(:,1),rec(:,2),'kx')
text(0.2,0.2,'K_\beta')
quiver(0.2,1-H/L1, 0.1,0,'MaxHeadsize', 30)
quiver(0,0,0.4*sin(th),0.4*cos(th),'MaxHeadsize', 30)


subplot(3,1,3)
%scatter(xx,zz,30,alpha,'filled'); hold on;
FV=griddata(xx,zz,alpha,X,Z);
imagesc(xa,za,FV); set(gca,'YDir','normal'); hold on;
plot(rec(:,1),rec(:,2),'kx')
text(0.2,0.2,'K_\alpha')
quiver(0.2,1-H/L1, 0.1,0,'MaxHeadsize', 30)
quiver(0,0,0.4*sin(th),0.4*cos(th),'MaxHeadsize', 30)

print('-djpeg',jpeg_file)
%file='angle_23/OUTPUT_FILES_KERNEL/proc000000_rhop_alpha_beta_kernel.dat';
%s=load(file);
%xx=s(:,1)/LENGTH;zz=s(:,2)/LENGTH;
%rhop=s(:,3); alpha=s(:,4); beta=s(:,5);
%field2=rhop;
% define surface grid
%maxv=max(abs(field))/3; ind=find(abs(field)>maxv); 
%field(ind)=maxv;
%xmax=max(xx); zmax=max(zz);
%xg = linspace(0,xmax,900);
%zg = linspace(0,zmax,300);
%[XI,YI] = meshgrid(xg,zg);
%FI = griddata(xx,zz,field,XI,YI);
%mesh(XI,YI,FI); view(2); 

%hold on;
% receivers


%subplot(3,1,2);
%scatter(xx,zz,30,field2+field1,'filled');hold on;

% plot assumed Moho at 40 km depth
%plot([0,3],[160,160]/L1,'k:')

% isochrons
% $$$ th=23*pi/180;
% $$$ beta=3.7; alpha=6.5;
% $$$ phi=linspace(0,pi);
% $$$ t0=5.0;
% $$$ r=t0./(1/beta-1/alpha*cos(th+pi/2-phi));
% $$$ r2=t0./(1/alpha-1/alpha*cos(th+pi/2-phi));
% $$$ plot(r.*cos(phi)/L1+rec(1),rec(2)-r.*sin(phi)/L1); hold on;
% $$$ plot(r2.*cos(phi)/L1+rec(1),rec(2)-r2.*sin(phi)/L1);
% $$$ axis([0,3,0,1])
