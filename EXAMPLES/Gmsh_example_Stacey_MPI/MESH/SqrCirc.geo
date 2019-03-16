// square
lc=0.1;
Point(1) = {1, 1, 0, lc};
Point(2) = {1, -1, 0, lc};
Point(3) = {-1, -1, 0, lc};
Point(4) = {-1, 1, 0, lc};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};

// big circle
R=0.2;
Centx=0.25;
Centy=0;
Point(5) = {Centx, R+Centy, 0, lc};
Point(6) = {R+Centx, Centy, 0, lc};
Point(7) = {-R+Centx, Centy, 0, lc};
Point(8) = {Centx, -R+Centy, 0, lc};
Point(9) = {Centx, Centy, 0, lc};
Circle(5) = {5, 9, 7};
Circle(6) = {7, 9, 8};
Circle(7) = {8, 9, 6};
Circle(8) = {6, 9, 5};
Line Loop(9) = {5, 6, 7, 8};
Line Loop(10) = {1, 2, 3, 4};

// small circle
R=0.1;
Centx=-.25;
Centy=0;
Point(105) = {Centx, R+Centy, 0, lc};
Point(106) = {R+Centx, Centy, 0, lc};
Point(107) = {-R+Centx, Centy, 0, lc};
Point(108) = {Centx, -R+Centy, 0, lc};
Point(109) = {Centx, Centy, 0, lc};
Circle(105) = {105, 109, 107};
Circle(106) = {107, 109, 108};
Circle(107) = {108, 109, 106};
Circle(108) = {106, 109, 105};
Line Loop(109) = {105, 106, 107, 108};

Plane Surface(10) = {10, 9, 109};
Plane Surface(11) = {9};
Plane Surface(12) = {109};

Recombine Surface{10,11,12};

// meshing
// quads mesh
Mesh.SubdivisionAlgorithm = 1;
Mesh.ElementOrder = 2;

//Mesh.Algorithm = 2; // frontal quad: 8;

// sizing
//Mesh.RecombinationAlgorithm = 2;
//Mesh.CharacteristicLengthFactor = 0.8;

Physical Line("Top") = {1};
Physical Line("Left") = {2};
Physical Line("Bottom") = {3};
Physical Line("Right") = {4};
Physical Surface("M1") = {10};
Physical Surface("M2") = {11,12};
