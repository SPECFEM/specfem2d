%
% function [bele1,bnod1,bele2,bnod2] = read_cubit_boundary(ifile,dnode,iele,nbound)
% Carl Tape, 19-Nov-2010
%
% INPUT
%   ifile   output file from cubit
%   dnode   nodes (vertices)
%   iele    described as a set of node indices
%
% OUTPUT
%   bele1   indices of boundary elements (use iele)
%   bnod1   indices into pairs of nodes describing each edge (use dnode)
%   bele2
%   bnod2
% 
% This assumes your model has TWO BOUNDARIES designated:
%   (1) free surface
%   (2) absorbing boundary
%
% CODE IN CUBIT TO GENERATE OUTPUT FILE (for Tromp2005 example):
% export Abaqus "Tromp2005_elements_nodes" overwrite cubitids
% block 1 curve 1
% block 2 curve 2 3 4
% export Abaqus "Tromp2005_boundaries" overwrite cubitids
%
% calls xxx
% called by xxx
%

function [bele1,bnod1,bele2,bnod2] = read_cubit_boundary(ifile,dnode,iele,nbound)

disp(sprintf('read_cubit_boundary.m: %i boundaries',nbound));
whos dnode iele

lines = textread(ifile,'%s','delimiter','\n');
n = length(lines);

nnod = length(dnode);
nele = length(iele);

jj = 0;
iedge = zeros(n,2);
%istart = []; kk = 0;
%iend = []; mm = 0;
istart = zeros(1,nbound); kk = 0;
for ii = 1:n
    ltemp = lines{ii};
    % if there are more than 8 characters in the line
    if length(ltemp) >= 8
        if strcmp(ltemp(1:8),'*ELEMENT')
            disp(ltemp);
            kk = kk+1;
            istart(kk) = ii;    % line number just before block
        else
            vals = sscanf(ltemp,'%f,');
            if length(vals)==3  % three columns of data (nodes or elements)
                %jj = jj+1;
                iedge(ii,:) = vals([2 3])';
            end
        end
    end
end

if isempty(istart), error('istart is empty'); end

% find starting and stopping parts for the two boundaries
istart1 = istart(1)+1;
if nbound==2
    istart2 = istart(2)+1;
    istop1 = istart2 - 2;
else
    istop1 = max(find(prod(iedge') ~= 0));
end
nele_b1 = istop1 - istart1 + 1;
bnod1 = iedge([istart1:istop1],:);
if any(bnod1(:)==0), error('zero entries in bnod1'); end

if nbound==2
    %istart2 = istart(2)+1;
    istop2 = max(find(prod(iedge') ~= 0));
    nele_b2 = istop2 - istart2 + 1;
    bnod2 = iedge([istart2:istop2],:);
    if any(bnod2(:)==0), error('zero entries in bnod2'); end
end

%------------------------------------------
% for each EDGE, determine the corresponding element
% NOTE: because these are boundary elements, there can only be ONE element
% corresponding to each edge.

bele1 = zeros(nele_b1,1);
for ii=1:nele_b1
    inds1 = bnod1(ii,:);
    [i1,~] = find(inds1(1) == iele);
    [i2,~] = find(inds1(2) == iele);
    bele1(ii) = intersect(i1,i2);
end

if nbound==2
    bele2 = zeros(nele_b2,1);
    for ii=1:nele_b2
        inds2 = bnod2(ii,:);
        [i1,~] = find(inds2(1) == iele);
        [i2,~] = find(inds2(2) == iele);
        %ii, inds2, i1, i2
        % NOTE (lqy): this explicitly requires that the boundary elements (edges) 
        % are on the boundary of the entire domain.
        bele2(ii) = intersect(i1,i2);
    end
end

ifig = 1;
if ifig==1
    figure; hold on;
    patch('faces',iele,'vertices',dnode,'facecolor','w');
    patch('faces',iele(bele1,:),'vertices',dnode,'facecolor','r');
    if nbound==2, patch('faces',iele(bele2,:),'vertices',dnode,'facecolor','c'); end
    title('boundary elements highlighted with color');
    axis equal, axis tight
end

%==========================================
