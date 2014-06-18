% read in forward seismograms

s=load('OUTPUT_FILES/S0005.AA.BXZ.semd');
t=s(:,1); disp=s(:,2);
tmin=12.0; tmax=21.0;
tshift=5.0; factor=0.2;
dt=t(2)-t(1);
ind=find(t>tmin & t<tmax);
ind_shift=floor(tshift/dt);

adj=zeros(size(disp));
adj(ind+ind_shift)=disp(ind)*factor;

plot(t,adj)

fid = fopen('S0005.AA.BXZ.adj','w');
fprintf(fid,'%6.2f %15.4g\n',[t';adj']); % the variable list has to be two rows and n columns
fclose(fid);