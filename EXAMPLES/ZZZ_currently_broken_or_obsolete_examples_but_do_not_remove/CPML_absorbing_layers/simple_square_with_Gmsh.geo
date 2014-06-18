//Square with Gmsh
lc=0.05;
//-----------------------------------------
p1=newp; Point(p1) = { 1,  1, 0, lc};
p2=newp; Point(p2) = { 1, -1, 0, lc};
p3=newp; Point(p3) = {-1, -1, 0, lc};
p4=newp; Point(p4) = {-1,  1, 0, lc};
//--------------------------------------------
l1=newl; Line(l1) = {p1, p4};
l2=newl; Line(l2) = {p4, p3};
l3=newl; Line(l3) = {p3, p2};
l4=newl; Line(l4) = {p2, p1};
//--------------------------------------------
ll0=newll; Line Loop(ll0) = {l1, l2, l3, l4};

s0=news; Plane Surface(s0) = {ll0};

Recombine Surface{s0};

Physical Line("Top") = {l1};
Physical Line("Left") = {l2};
Physical Line("Bottom") = {l3};
Physical Line("Right") = {l4};
Physical Surface("M1") = {s0};
