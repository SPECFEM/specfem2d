/* This routine to be called from fortran program writes
 * the ParaView image file (.vti) 
 * Hom Nath Gharti*/
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

int FC_FUNC_(write_vtk_imagec,WRITE_VTK_IMAGEC)(
  int *ntframe_f,int *nx_f,int *nz_f,float *dx_f,float *dz_f,int *it_f,float *dt_f, 
  int *nstep_f,float data[*nz_f][*nx_f],char *vtag)
{
	int i,j,nvar;
	int nx,ny,nz,it,nstep,ntframe;
	int off[10];
	int bytes[10];
	float dt,dx,dy,dz,ox,oy,oz,temp;
	char fname[80],outfname[80];
	FILE *outf0,*outf1;

	nx=*nx_f;
	nz=*nz_f;
	it=*it_f;
	
	dt=*dt_f;
	nstep=*nstep_f;
  ntframe=*ntframe_f;
	ox=oy=oz=0.0;
	dx=*dx_f; dz=*dz_f;
	dy=dx; /* For ZX plane it shouldn't matter */
	nvar=1; // Plot only one variable set in each file
	off[0]=0;
	for (i=0; i<nvar; i++){
		bytes[i]=(nx*nz)*sizeof(float);
		if(i<nvar-1)off[i+1]=off[i]+sizeof(int)+bytes[i];
		bytes[i]=bytes[i]+sizeof(int);
	}

  if(it != -9999 && nstep != -9999){
    sprintf(outfname,"./OUTPUT_FILES/%s.pvd",vtag);
    if(it==1 || it==ntframe){
      outf0=fopen(outfname,"w");
      fprintf(outf0,"<?xml version=\"1.0\"?>\n");
      fprintf(outf0,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
      fprintf(outf0,"  <Collection>\n");	
    }else{
      outf0=fopen(outfname,"a");
    }
  }
		
  if(it != -9999 && nstep != -9999){
	  sprintf(outfname,"./OUTPUT_FILES/%s_%07d.vti",vtag,it);
	  sprintf(fname,"%s_%07d.vti",vtag,it);
	  fprintf(outf0,"    <DataSet timestep=\"%f\" part=\"001\" file=\"%s\" name=\"Asmb:FRAME\"/>\n",it*dt,fname);
  }else{
	  sprintf(outfname,"./OUTPUT_FILES/%s.vti",vtag);
  }
	/* printf("%s\n",outfname); */	
	outf1=fopen(outfname,"wb");

	fprintf(outf1,"<?xml version=\"1.0\"?>\n");
	fprintf(outf1,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(outf1,"  <ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%f %f %f\" Spacing=\"%f %f %f\">\n",0,nx-1,0,0,0,nz-1
		,ox,oy,oz
		,dx,dy,dz);
	fprintf(outf1,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",0,nx-1,0,0,0,nz-1);
	fprintf(outf1,"      <PointData Scalars=\"%s\">\n",vtag);
	fprintf(outf1,"        <DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"%d\" />\n",vtag,off[0]);

	fprintf(outf1,"      </PointData>\n");
	fprintf(outf1,"      <CellData>\n");
	fprintf(outf1,"      </CellData>\n");
	fprintf(outf1,"    </Piece>\n");
	fprintf(outf1,"  </ImageData>\n");
	fprintf(outf1,"  <AppendedData encoding=\"raw\">\n");
	fprintf(outf1,"_");
		
	fwrite(&bytes[0],sizeof(int),1,outf1);
	/* fwrite(data,sizeof(data),1,outf1); *//* Need to check for this short statement, why is it not working? */ 
	for(j=0; j<nz; j++){		
		for(i=0; i<nx; i++){
			temp=data[j][i];
			fwrite(&temp,sizeof(float),1,outf1);
		}
	}
	fprintf(outf1,"\n");
	fprintf(outf1,"  </AppendedData>\n");
	fprintf(outf1,"</VTKFile>\n");
	fclose(outf1);

  if(it != -9999 && nstep != -9999){
    if(it==nstep){
      fprintf(outf0,"  </Collection>\n");
      fprintf(outf0,"</VTKFile>\n");
    }
	  fclose(outf0);
  }
  return(1);
}
