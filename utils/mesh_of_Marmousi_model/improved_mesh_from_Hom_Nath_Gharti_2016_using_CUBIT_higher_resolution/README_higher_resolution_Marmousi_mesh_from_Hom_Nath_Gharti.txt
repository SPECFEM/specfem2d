Subject: Re: Marmousi / here is the mesh created by Yann with Gmsh
Date: Thu, 12 May 2016 12:27:32 -0400
From: Hom Nath Gharti
To: Dimitri Komatitsch
CC: Yann Capdeville, Bence Solymosi, Jeannot Trampert, Jeroen Tromp, Ryan T. Modrak, Tarje Nissen-Meyer, Vadim Monteiller, Yang Luo

Hi Dimitri,

Sorry for delay! I finally had time to work on this.

Attached complete mesh files and TRELIS/CUBIT files. Mesh quality is
now improved from scaled jacobian 1.083e-02 to
 6.063e-02.  Those files are ready to run with specfem2d but still
missing are Par_file and SOURCE files. As Yann previously mentioned,
material properties file needs to be decoded! I did not delete or
merge any surfaces, therefore the original material properties file
should work as is once decoded to specfem2d format.

I put the origin on the top left corner, define the absorbing
boundaries on left, right and bottom, and free surface on top.

Please note that the .trelis and .cub files are not compatible with
old Cubit version. Since smoothing/blending/chamferring sharp corners
is a manual process, creating a Journal file is unrealistic.

Let me know if you find any issue/s.

Best,
Hom Nath

