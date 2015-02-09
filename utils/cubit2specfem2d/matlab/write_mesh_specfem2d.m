%
% function write_mesh_specfem2d
% Carl Tape, 07-Jan-2010
%
% This writes 5 or 6 mesh files that are read into SPECFEM2D.
%
% NOTE: SPECFEM2D expects ordering COUNTER-clockwise, but CUBIT may output
%       either convention, even within the same mesh (WHY?). This function
%       ensures that the output is consistently COUNTER-clockwise.
% 
% calls xxx
% called by xxx
%

function write_mesh_specfem2d(otag,dnode,iele,mele,bele1,bnod1,bele2,bnod2)

disp('entering write_mesh_specfem2d.m');

% surface 1: free (may be empty)
ifree = 1;
iabs = 1;
if and(isempty(bele1),isempty(bnod1))
    disp('no free surface boundary');
else
    disp('free surface boundary present');
    nele_b1 = length(bele1);
    % surface 1: free
    ofile = [otag '_surface_free'];
    disp(['writing ' ofile]);
    fid = fopen(ofile,'w');
    fprintf(fid,'%i\n',nele_b1); 
    for ii = 1:nele_b1
        fprintf(fid,'%14i%14i%14i%14i\n',bele1(ii),2,bnod1(ii,1),bnod1(ii,2));   
    end
    fclose(fid);
end

% surface 2: absorbing (may be empty)
if and(isempty(bele2),isempty(bnod2))
    disp('no absorbing boundary');
else
    disp('absorbing boundary present');
    nele_b2 = length(bele2);
    ofile = [otag '_surface_absorb'];
    disp(['writing ' ofile]);
    fid = fopen(ofile,'w');
    fprintf(fid,'%i\n',nele_b2); 
    for ii = 1:nele_b2
        fprintf(fid,'%14i%14i%14i%14i\n',bele2(ii),2,bnod2(ii,1),bnod2(ii,2));   
    end
    fclose(fid);
end

nnod = length(dnode);
nele = length(iele);

% nodes
ofile = [otag '_nodes'];
disp(['writing ' ofile]);
fid = fopen(ofile,'w');
fprintf(fid,'%16i\n',nnod); 
for ii = 1:nnod
    fprintf(fid,'%18.10e%18.10e\n',dnode(ii,1),dnode(ii,2));   
end
fclose(fid);

% mesh of elements
% NOTE: SPECFEM2D expects ordering COUNTER-clockwise, but CUBIT may output
%       either convention, even within the same mesh (WHY?).
ofile = [otag '_elements'];
disp(['writing ' ofile]);
fid = fopen(ofile,'w');
fprintf(fid,'%i  4\n',nele); 
for ii = 1:nele
    iv = iele(ii,:);
    if iv > nnod
        nnod, iv, ii, iele(ii,:)
        error('invalid index into nodes list');
    end
    x1 = dnode(iv,1);
    y1 = dnode(iv,2);
    if ispolycw(x1,y1)
        iv = iele(ii,[1 4 3 2]);
    end
    fprintf(fid,'%i %i %i %i\n',iv);
end
fclose(fid);

% material index for each element
ofile = [otag '_material'];
disp(['writing ' ofile]);
fid = fopen(ofile,'w');
for ii = 1:nele
    fprintf(fid,'%i\n',mele(ii));   
end
fclose(fid);

%==========================================
