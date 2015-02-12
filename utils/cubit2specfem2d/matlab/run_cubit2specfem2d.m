%
% test_specfem2d.m
% Carl Tape, 01-June-2010
%
% 
%
%

clc
clear
close all
format short, format compact

iwrite = 1;

%=========================================================================
% read cubit abaqus file and generate text files for SPECFEM2D

sfac = 1;
icol = [2 4];
nbound = 2;
nnod = 2673; nele = 2560; stag = 'Tromp2005'; iunit = [1];
dir0 = '/data/svn/seismo/2D/SPECFEM2D/UTILS/cubit2specfem2d/matlab/';

ifile1 = [dir0 stag '_elements_nodes'];
ifile2 = [dir0 stag '_boundaries'];
[dnode,iele,mele] = read_cubit_basic(ifile1,nele,nnod,iunit,icol);
[bele,bnod,bele2,bnod2] = read_cubit_boundary(ifile2,dnode,iele,nbound);

if iwrite==1
    otag = [dir0 stag];
    write_mesh_specfem2d(otag,dnode*sfac,iele,mele,bele,bnod,bele2,bnod2);
end

%==========================================================================
