close all; clear all;
% input parameters
H=40; Vp1=6.5; Vp2=8.3; Vs1=3.7; Vs2=4.7; X0=-150;
dir='OUTPUT_FILES_FORWARD/';

% receivers
Xs=100;Xe=500;
nsta=30;
if nsta == 1
  xarray=[Xs];
else
  xarray=linspace(Xs,Xe,nsta);
end

% incident waves
angle_inc=18;
th=angle_inc*pi/180; p=sin(th)/Vp2;

tmin=zeros(nsta,1); tmax=zeros(nsta,1);
[tP, tPs, tPpp, tPps, tPss]=ttplane(X0,xarray,p,H,Vp2,Vs2,Vp1,Vs1);

% reduced time
rt=tP;
% plot distance between traces
dy=0.5;
% start of seismograms

%%% select t1 and t2 of adjoint source
t1=2; t2=7; rtaper=0;  % select Ps phase
t1=2; t2=30; 

comp = ['BXZ'; 'BXX'];
% plot seismograms
for i = 1 : nsta
  X=(Xe-Xs)/(nsta-1)*(i-1)+Xs;
  
  for j = 1 : 2
    figure(j);
    file=sprintf([dir 'S00%02d.AA.' comp(j,:) '.semd'],i);
    s=load(file);
    ym = (i-1)*dy;
    plot(s(:,1)-rt(i),s(:,2)+ym); hold on; grid on;
    
    % pick start and end of adjoint seismograms (after aligning P)
    tmin(i)=t1; 
    tmax(i)=min(t2,s(end,1)-rt(i)); 
    plot(tmin(i),ym,'rx',tmax(i),ym,'mx');
    
    % generate adjint source
    ind=find(s(:,1)-rt(i)>tmin(i) & s(:,1)-rt(i)<tmax(i));
    sadj=zeros(size(s(:,2)));
    taper=tukeywin(length(ind),rtaper);
    sadj(ind)=s(ind,2).*taper;
    plot(s(ind,1)-rt(i),sadj(ind)+ym,'g');
    
    % write adjoint source
    fid = fopen(strcat(file,'.adj'),'w');
    fprintf(fid,'%6.2f %12.4g\n',[s(:,1)';sadj']);
    fclose(fid);
  end
  
  file=sprintf([dir 'S00%02d.AA.BXY.semd'],i);
  fid = fopen(strcat(file,'.adj'),'w');
  fprintf(fid,'%6.2f %12.4g\n',[s(:,1)';zeros(size(s(:,1)))']);
  fclose(fid);
end

for j = 1 : 2
  figure(j)
  plot(tPs-tP,(0:nsta-1)*dy); text(tPs(1)-tP(1),-0.5,'tPs');
  plot(tPss-tP,(0:nsta-1)*dy); text(tPss(1)-tP(1),-0.5,'tPss');
  plot(tPps-tP,(0:nsta-1)*dy); text(tPps(1)-tP(1),-0.5,'tPps');
  plot(tPpp-tP,(0:nsta-1)*dy); text(tPpp(1)-tP(1),-0.5,'tPpp');
end
