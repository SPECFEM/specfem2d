
%% load all curves
S1=load('baty_my_salt.txt');
S2=load('c1_my_salt.txt');
S3=load('c2_my_salt.txt');
S4=load('c3_my_salt.txt');
S5=load('c4_my_salt.txt');
S6=load('c5_my_salt.txt');
S7=load('c6_my_salt.txt');
S8=load('dome_my_salt.txt');


%% define model zone

S8_bis=flipud(S8);
% define different zones (close curves) 

% ocean 
zone0=[S1;200 0; 0 0; S1(1,:)];

% layer 1
zone1=[S1;flipud(S2);S1(1,:)];

%layer 2
zone2=[S2; flipud(S3);S2(1,:)];

% layer 3
zone3=[S3; flipud(S7);S8_bis(10:1:29,:);flipud(S6);S3(1,:)];

% layer 4
zone4=[S6;S8_bis(30:1:33,:);flipud(S4);0 70;S6(1,:)];

%layer 5
zone5=[S4;S8_bis(34:1:end,:);S4(1,:)];

% layer 6 
zone6=[S7;200 70;flipud(S5);S8_bis(6:1:9,:);S7(1,:)];

% layer 7
zone7=[S8_bis(1:1:5,:);S5;S8_bis(1,:)];

% salt dome
zone8=[S8; S8(1,:)];

figure; hold on
plot(zone1(:,1),-zone1(:,2),'.-')
plot(zone2(:,1),-zone2(:,2),'.-k')
plot(zone3(:,1),-zone3(:,2),'.-r')
plot(zone4(:,1),-zone4(:,2),'.-k')
plot(zone5(:,1),-zone5(:,2),'.-r')
plot(zone6(:,1),-zone6(:,2),'.-')
plot(zone7(:,1),-zone7(:,2),'.-r')
plot(zone8(:,1),-zone8(:,2),'.-y')


%% grid model 
xout=-100000;
dh=0.1;
x=0:dh:200;
z=0:dh:70;
nz=length(z);
nx=length(x);
M=zeros(nx,nz);
for iz=1:1:nz;
    for ix=1:1:nx
        X=x(ix);
        Z=z(iz);
        
         [ in ] = Inside_zone(X,Z,zone0,xout);
         if (mod(in,2)==1); M(ix,iz)=1;end;
         
         [ in ] = Inside_zone(X,Z,zone1,xout);
         if (mod(in,2)==1); M(ix,iz)=2;end;
         
         [ in ] = Inside_zone(X,Z,zone2,xout);
          if (mod(in,2)==1); M(ix,iz)=3;end;
         
         [ in ] = Inside_zone(X,Z,zone3,xout);
          if (mod(in,2)==1); M(ix,iz)=4;end;
         
         [ in ] = Inside_zone(X,Z,zone4,xout);
          if (mod(in,2)==1); M(ix,iz)=5;end;
         
         [ in ] = Inside_zone(X,Z,zone5,xout);
          if (mod(in,2)==1); M(ix,iz)=6;end;
         
         [ in ] = Inside_zone(X,Z,zone6,xout);
          if (mod(in,2)==1); M(ix,iz)=7;end;
         
          
         [ in ] = Inside_zone(X,Z,zone7,xout);
          if (mod(in,2)==1); M(ix,iz)=8;end;
         
          
         [ in ] = Inside_zone(X,Z,zone8,xout);
          if (mod(in,2)==1); M(ix,iz)=9;end;
         
    end
end
    
% remove remaining 0 

for iz=2:1:nz-1
    for ix=2:1:nx-1
        if (M(ix,iz)==0);
            found=0;
            for jj=-1:1:1
                for ii=-1:1:1
                    if (M(ix+ii,iz+jj) >0); M(ix,iz)= M(ix+ii,iz+jj);found=1;break ;end;
                end 
                if (found); break;end
            end 
           
        end
    end
end

% remove zeors on boundary
M(1,:)=M(2,:);
M(:,1)=M(:,2);
M(nx,:)=M(nx-1,:);
M(:,nz)=M(:,nz-1);

figure;
imagesc(M');axis image;


%% prolongement de M
p=20;
Mp=zeros(nx+2*p,nz+2*p);
Mp(p+1:1:nx+p,p+1:1:nz+p)=M;

for ip=1:1:p+1
    Mp(ip,:)=Mp(p+1,:);
    Mp(:,ip)=Mp(:,p+1);
end

for ip=1:1:p+1
    Mp(nx+p+ip-1,:)=Mp(nx+p,:);
    Mp(:,nz+p+ip-1)=Mp(:,nz+p);
end


%% define physical properties 

Vp=[1500 2000 2500 2800 3000 3500 3000 3500 4700];
Vs=[0 1176 1470 1555 1875 2100 1875 2100 2880];
rho=[1000 2000 2500 2800 3000 3500 3000 3500 3900 ];
Qkappa=[100 100 100 100 100 100 100 100 100];
Qmu=[100 100 100 100 100 100 100 100 100];

% 1  1 1000.d0 1500.d0    0.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 2  1 2000.d0 2000.d0 1176.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 3  1 2500.d0 2500.d0 1470.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 4  1 2800.d0 2800.d0 1555.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 5  1 3000.d0 3000.d0 1875.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 6  1 3500.d0 3500.d0 2100.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 7  1 3900.d0 4700.d0 2880.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 8  1 1000.d0 1500.d0    0.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0
% 9  1 2000.d0 2000.d0 1176.0d0  0 0 100.d0 100.d0 0 0 0 0 0 0


nnx=size(Mp,1);
nnz=size(Mp,2);
for iz=1:1:nnz
    for ix=1:1:nnx
        Mvp(ix,iz)=Vp(Mp(ix,iz));
        Mvs(ix,iz)=Vs(Mp(ix,iz));
        Mrho(ix,iz)=rho(Mp(ix,iz));
    end 
end 

