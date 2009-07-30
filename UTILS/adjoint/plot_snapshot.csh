gmtset BASEMAP_TYPE plain ANOT_FONT_SIZE 13 HEADER_FONT_SIZE 8
makecpt -Cseis -T-1.110e-07/1.110e-07/2.018e-08 > color.cpt
sed 's/^N.*/N       255     255     255/' color.cpt > color1.cpt
sed 's/^B.*/B       255     255     255/' color1.cpt > color.cpt
echo OUTPUT_FILES/snapshot_rho_kappa_mu_0003000
psxy -JX6i/2i -R0/1/0/1 -X3 -Y-3 -K -V -P <<EOF >plot_01.ps
EOF
awk '{print $1,$2,$3 * 1}' OUTPUT_FILES/snapshot_rho_kappa_mu_0003000 | pscontour -JX6i/2i -R0/1/0/1 -B0.5/0.1:."rho":WeSn -A- -Ccolor.cpt -I -K -O -V -Y8 -P>> plot_01.ps
echo Fig ok
psscale -Ccolor.cpt -D3i/-0.5i/5c/0.1h -B9.99e-06 -K -O -P >> plot_01.ps 
awk '{print $1,$2,$4 * 1}' OUTPUT_FILES/snapshot_rho_kappa_mu_0003000 | pscontour -JX6i/2i -R0/1/0/1 -B0.5/0.1:."kappa":WeSn -A- -Ccolor.cpt -I -K -O -V -Y8 -P>> plot_01.ps
echo Fig ok
awk '{print $1,$2,$5 * 1}' OUTPUT_FILES/snapshot_rho_kappa_mu_0003000 | pscontour -JX6i/2i -R0/1/0/1 -B0.5/0.1:."mu":WeSn -A- -Ccolor.cpt -I -K -O -V -Y8 -P>> plot_01.ps
echo Fig ok
pstext -JX -R -K -O -P -N >>plot_01.ps<<EOF
 0.6 0.45 20 0 4 RM time = ? s
EOF
psxy -JX -R -O -P -V <<EOF>>plot_01.ps
EOF
