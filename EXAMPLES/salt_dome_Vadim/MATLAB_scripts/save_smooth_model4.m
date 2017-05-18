


nnx=size(Mevp,1);
nnz=size(Mevp,2);


xmin=-ddhx*dh;
zmin=0; %(nnz-1)*dh;

x=xmin:dh:(nnx-1)*dh + xmin;
z=zmin:dh:(nnz-1)*dh + zmin;

figure;
imagesc(x/10,(z)/10,Mevp');axis image;

A_to_save=zeros(nnz*nnx,7);
k=0;
for iz=nnz:-1:1
    nnz-iz+1 
    for ix=1:1:nnx
        k=k+1;
        A_to_save(k,:)= [x(ix)*100 -z(iz)*100 Mevp(ix,iz) Mevs(ix,iz)  Merho(ix,iz) 100 100 ];
    end
end

vpmin=min(min(Mevp));
vpmax=max(max(Mevp));
vsmin=min(min(Mevs));
vsmax=max(max(Mevs));
rhomin=min(min(Merho));
rhomax=max(max(Merho));

fid=fopen('model_init_smooth_200m','w');
fprintf(fid,' %f %f %f %f \n',[min(x*100) max(x*100) min(-z*100) max(-z*100)]);
fprintf(fid,' %f %f \n',[dh*100 dh*100]);
fprintf(fid,' %d %d \n',[nnx nnz]);
fprintf(fid,' %f %f %f %f %f %f \n',[vpmin vpmax vsmin vsmax rhomin rhomax]);
fprintf(fid,'%f %f %f %f %f %f %f\n',A_to_save')

fclose(fid)