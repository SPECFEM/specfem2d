%
% function [dnode,iele,mele] = read_cubit_basic(ifile,nele,nnod,iunit,icol1)
% Carl Tape, 19-Nov-2010
% 
% This function reads an abaqus text file output from cubit to obtain the
% nodes and elements.
%
% INPUT
%   ifile   input file
%   nele    number of elements (listed in cubit window)
%   nnod    number of nodes (listed in cubit window)
%   iunit   index vector for blocks (NOTE: order must match what is
%               assigned in the SPECFEM2D Par_file)
%   icol1   OPTIONAL: denotes 2D mesh, indices pick the columns you want
%
% OUTPUT
%   dnode   nnod x 2 matrix of node vertices
%   iele    nele x 4 matrix of indices into nodes
%   mele    nele x 1 vector of material properties
%
% calls xxx
% called by test_specfem2d.m
%

function [dnode,iele,mele] = read_cubit_basic(ifile,nele,nnod,iunit,icol1)

nunit = length(iunit);

if nargin==5    % 2D mesh
    ncol1 = 2; ncol2 = 4; i2D = 1;
else            % 3D mesh
    ncol1 = 3; ncol2 = 8; icol1 = [2 3 4]; i2D = 0;
end
icol2 = 2:(ncol2+1);
disp(sprintf('read_cubit_basic.m: %s', ifile));

lines = textread(ifile,'%s','delimiter','\n');

% line numbers for start and finish
in1 = 10;
in2 = in1 + nnod - 1;
ie1 = in2 + 3;
ie2 = ie1 + nele + nunit - 1;

% read in nodes
dnode = zeros(nnod,ncol1);
inode = zeros(nnod,1);
for ii = 1:nnod
    kk = in1-1+ii;
    vals = sscanf(lines{kk},'%f,');
    dnode(ii,:) = vals(icol1)';
    inode(ii) = vals(1);     % NOTE: indexing may not run from 1:nnod
end

iindex = 0;
disp(sprintf('%i nodes, imin = %i, imax = %i',nnod,min(inode),max(inode)));
if max(inode) ~= nnod
    iindex = 1;
    disp('WARNING: node indexing from cubit does not run from 1 to nnod');
    disp('try using the cubit command: compress node all');
    error('exit');
end

% read in elements
iele = zeros(nele,ncol2);
mele = zeros(nele,1); im = 0;
ieleind = zeros(nele,1);  % index listed by cubit
jj = 0;
imax = ie2-ie1+1;
for ii = 1:imax
    kk = ie1-1+ii;
    ltemp = lines{kk};
    if strcmp(ltemp(1),'*')
        disp(ltemp);
        im = im+1;  % advance to next block of elements
    else
        vals = sscanf(ltemp,'%f,');
        jj = jj+1;
        ieleind(jj) = vals(1);
        iele(jj,:) = vals(icol2)';
        
        % assign material index
        % WARNING: the order of blocks exported by cubit does not
        % necessarily correspond to the index order of the surfaces
        %mele(jj) = im;
        mele(jj) = iunit(im);
    end
end
if any(iele(:)==0), error('zero entries in iele'); end

whos dnode iele mele

ifig = 1;
if ifig==1
    if i2D==1       % 2D mesh
        figure;
        patch('faces',iele,'vertices',dnode,'facecolor','flat','FaceVertexCData',[1:nele]');
        colorbar; title('colored in order of element index');
        axis equal, axis tight

        figure;
        patch('faces',iele,'vertices',dnode,'facecolor','flat','FaceVertexCData',mele);
        colorbar; title('colored by material index');
        axis equal, axis tight
        
    else            % 3D mesh
        % all nodes
        figure;
        plot3(dnode(:,1), dnode(:,2), dnode(:,3), '.');
        axis equal; box on;
        
        % plot a particular element
        ipick = 1;
        pele = iele(ipick,:);
        figure;
        plot3(dnode(pele,1), dnode(pele,2), dnode(pele,3), 'ro','markersize',12,'markerfacecolor','w');
        for ii=1:8
            stlab = sprintf('(%i, %i)',ii,pele(ii));
            text(dnode(pele(ii),1), dnode(pele(ii),2), dnode(pele(ii),3),stlab);
        end
        axis equal; box on;
        title(sprintf('%i nodes defining element %i (cubit index %i)',ncol2,ipick,ieleind(ipick)));
    end
end

%==========================================
