close all; clear all
sframe=200; eframe=5000; dframe=sframe
i=0;
for iframe = [5, sframe:dframe:eframe]
  i = i+1
  file=sprintf('OUTPUT_FILES/wavefield%07d_02_000.txt',iframe);
  s=load(file);
  if (i == 1)
    x=s(:,1); z=s(:,2);
  end
  sx=s(:,3); sz=s(:,5);
  figure(i)
  scatter(x,z,30,sz);  % or sx
end