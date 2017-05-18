

%% prolongement de M
p=500;
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
Vs=[1176 1176 1470 1555 1875 2100 1875 2100 2880];
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

%% smooth model 

% gaussian kernel
s=20;n=10*s;
im=n/2;
G=zeros(n,n);
sx2=s*s;
sy2=sx2;

for iy=1:1:n
    for ix=1:1:n
        G(ix,iy) = exp( -0.5*( (ix-im)^2/sx2 + (iy-im)^2/sy2 ) ) / (2*pi*sx2);
    end 
end

ceof=sum(sum(G));

Ms=conv2(Mp,G,'same')/ceof;
Msvp=conv2(Mvp,G,'same')/ceof;
Msvs=conv2(Mvs,G,'same')/ceof;
Msrho=conv2(Mrho,G,'same')/ceof;

% extract model 
ddhz=100;
ddhx=200;
Me=Ms(s+p-ddhx:1:nx+s+p-1+ddhx, s+p:1:ddhz+nz+s+p-1);
% figure;
% imagesc(Me');axis image;

Mevp=Msvp(p-ddhx:1:nx+p-1+ddhx, p:1:ddhz+nz+p-1);
Mevs=Msvs(p-ddhx:1:nx+p-1+ddhx, p:1:ddhz+nz+p-1);
Merho=Msrho(p-ddhx:1:nx+p-1+ddhx, p:1:ddhz+nz+p-1);

Mvp1=Mvp(p-ddhx:1:nx+p-1+ddhx, p:1:ddhz+nz+p-1);

figure;
imagesc(Mevp');axis image;
figure;
imagesc(Mvp1');axis image;



