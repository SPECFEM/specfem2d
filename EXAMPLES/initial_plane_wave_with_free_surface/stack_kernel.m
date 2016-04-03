close all; clear all;
% choose incident angles and directions to stack
LENGTH=250000; L1=LENGTH/1000; H=40;
nrec=30;
% modify directions of incident waves for stacked kernels
dir='+-'; ndir=length(dir);
angle=[18,20,23];

for idir = 1 : ndir
  cdir=dir(idir);
  for j = 1 : length(angle)
    cang=angle(j);
    
    kdir=sprintf('ang-%2.2i/nrec-%2.2i-dir-%s/',cang,nrec,cdir)
    file=[kdir 'OUTPUT_FILES_KERNEL/proc000000_rhop_alpha_beta_kernel.dat'];
    s=load(file);
    if idir == 1 & j == 1
      rhop=s(:,3); alpha=s(:,4); beta=s(:,5);
      file=[kdir 'OUTPUT_FILES_FORWARD/receiver.dat']; rec=load(file);
      rec=rec/LENGTH; nrec=length(rec);
    else
      rhop=rhop+s(:,3); alpha=alpha+s(:,4); beta=beta+s(:,5);
    end
  end
end

xx=s(:,1)/LENGTH;zz=s(:,2)/LENGTH;
minx=min(xx); maxx=max(xx); minz=min(zz); maxz=max(zz);
nx=400; nz=400;
% define mesh
xa=linspace(minx,maxx,nx);za=linspace(minz,maxz,nz);
[X,Z] = meshgrid(xa,za);

subplot(3,1,1)
FV=griddata(xx,zz,rhop,X,Z);
imagesc(xa,za,FV); set(gca,'YDir','normal'); hold on;
plot(rec(:,1),rec(:,2),'kx')
text(0.2,0.2,'K_{\rho_p} ')
quiver(0.2,1-H/L1, 0.1,0,'MaxHeadsize', 30)
title(['stacked by inc angle & dir=',num2str(ndir), ', nrec=', num2str(nrec)]);

subplot(3,1,2)
FV=griddata(xx,zz,beta,X,Z);
imagesc(xa,za,FV); set(gca,'YDir','normal'); hold on;
plot(rec(:,1),rec(:,2),'kx')
plot(rec(:,1),rec(:,2),'kx')
text(0.2,0.2,'K_\beta')
quiver(0.2,1-H/L1, 0.1,0,'MaxHeadsize', 30)

subplot(3,1,3)
FV=griddata(xx,zz,alpha,X,Z);
imagesc(xa,za,FV); set(gca,'YDir','normal'); hold on;
plot(rec(:,1),rec(:,2),'kx')
text(0.2,0.2,'K_\alpha')
quiver(0.2,1-H/L1, 0.1,0,'MaxHeadsize', 30)

jpeg=['stacked-kernel-' num2str(nrec) '-' num2str(ndir) '.jpg']
print('-djpeg',jpeg)