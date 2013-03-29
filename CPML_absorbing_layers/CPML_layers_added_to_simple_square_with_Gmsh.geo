// PML addition to a homogeneous square

Include "simple_square_with_Gmsh.geo";

// Define the PML structure
// Only two configurations are implemented at the moment
// 1 : PML everywhere
// 2 : Only Top free 

ConfigPML=1;

NbElPML = 3;
NbTrans = 4;
//---------------------------------------------------------------------------------------------------
p11=newp; Translate {NbElPML*lc, 0, 0}            { Duplicata{ Point{p1}; } }
p12=newp; Translate {NbElPML*lc, NbElPML*lc, 0}   { Duplicata{ Point{p1}; } }
p13=newp; Translate {0, NbElPML*lc, 0}            { Duplicata{ Point{p1}; } }
//
p21=newp; Translate {NbElPML*lc, 0, 0}            { Duplicata{ Point{p2}; } }
p22=newp; Translate {NbElPML*lc, -NbElPML*lc, 0}  { Duplicata{ Point{p2}; } }
p23=newp; Translate {0, -NbElPML*lc, 0}           { Duplicata{ Point{p2}; } }
//
p33=newp; Translate {-NbElPML*lc, 0, 0}           { Duplicata{ Point{p3}; } }
p32=newp; Translate {-NbElPML*lc, -NbElPML*lc, 0} { Duplicata{ Point{p3}; } }
p31=newp; Translate {0, -NbElPML*lc, 0}           { Duplicata{ Point{p3}; } }
//
p41=newp; Translate {-NbElPML*lc, 0, 0}           { Duplicata{ Point{p4}; } }
p42=newp; Translate {-NbElPML*lc, NbElPML*lc, 0}  { Duplicata{ Point{p4}; } }
p43=newp; Translate {0, NbElPML*lc, 0}            { Duplicata{ Point{p4}; } }

If (ConfigPML==2)

	// Configuration Top Free

	//------------Right---------------
	lr1=newl; Line(lr1) = {p2,p21};
	lr2=newl; Line(lr2) = {p21,p11};
	lr3=newl; Line(lr3) = {p11,p1};

	Transfinite Line{lr1,lr3} = NbTrans;

	ll1=newll; Line Loop(ll1) = {lr1,lr2,lr3,-l4};
	s1=news; Plane Surface(s1) = {ll1};
	Transfinite Surface{s1};
	Physical Surface("R") = {s1};

	//----------Bottom Right----------------
	//l1=newl; Line(l1) = {2,p2};
	lr4=newl; Line(lr4) = {p2,p23};
	lr5=newl; Line(lr5) = {p23,p22};
	lr6=newl; Line(lr6) = {p22,p21};

	Transfinite Line{lr4,lr5,lr6} = NbTrans;

	ll2=newll; Line Loop(ll2) = {-lr1,lr4,lr5,lr6};
	s2=news; Plane Surface(s2) = {ll2};
	Transfinite Surface{s2};
	Physical Surface("RB") = {s2};
	//Transfinite Surface{s1} = {1,p1,p2,2};

	//---------Bottom----------------------
	//l1=newl; Line(l1) = {2,p3};
	lr7=newl; Line(lr7) = {p3,p31};
	lr8=newl; Line(lr8) = {p31,p23};

	Transfinite Line{lr7} = NbTrans;

	ll3=newll; Line Loop(ll3) = {-lr4,-l3,lr7,lr8};
	s3=news; Plane Surface(s3) = {ll3};
	Transfinite Surface{s3};
	Physical Surface("B") = {s3};

	//-----Bottom Left-----------------
	//l1=newl; Line(l1) = {3,p5};
	lr9=newl;  Line(lr9)  = {p3,p33};
	lr10=newl; Line(lr10) = {p33,p32};
	lr11=newl; Line(lr11) = {p32,p31};

	Transfinite Line{lr9,lr10,lr11}=NbTrans;

	ll4=newll; Line Loop(ll4) = {-lr7,lr9,lr10,lr11};
	s4=news; Plane Surface(s4) = {ll4};
	Transfinite Surface{s4};
	Physical Surface("LB") = {s4};

	// -------------Left---------------------------
	//l1=newl; Line(l1) = {3,p6};
	lr12=newl; Line(lr12) = {p4,p41};
	lr13=newl; Line(lr13) = {p41,p33};

	Transfinite Line{lr12}=NbTrans;

	ll5=newll; Line Loop(ll5) = {-lr9,-l2, lr12,lr13};
	s5=news; Plane Surface(s5) = {ll5};
	Transfinite Surface{s5};
	Physical Surface("L") = {s5};
	//
	//Transfinite Surface "*";
	//Physical Surface("PML") = {s1,s2,s3,s4,s5};
	// Bords du domaine PML absorbants
	Physical Line ("Left_PML") = {lr10,lr13};
	Physical Line ("Right_PML") = {lr2,lr6};
	Physical Line ("Bottom_PML") = {lr5,lr8,lr11};
	// Ici Bord du domaine PML free continuité du domaine initial
	Physical Line ("Top_PML") = {lr3,lr12,l1};  //+Top Free
EndIf

If (ConfigPML==1)

	// Configuration Aborbant partout

	//------------Right---------------
	lr1=newl; Line(lr1) = {p2,p21};
	lr2=newl; Line(lr2) = {p21,p11};
	lr3=newl; Line(lr3) = {p11,p1};

	Transfinite Line{lr1,lr3} = NbTrans;

	ll1=newll; Line Loop(ll1) = {lr1,lr2,lr3,-l4};
	s1=news; Plane Surface(s1) = {ll1};
	Transfinite Surface{s1};
	Physical Surface("R") = {s1};

	//----------Bottom Right----------------
	lr4=newl; Line(lr4) = {p2,p23};
	lr5=newl; Line(lr5) = {p23,p22};
	lr6=newl; Line(lr6) = {p22,p21};

	Transfinite Line{lr4,lr5,lr6} = NbTrans;

	ll2=newll; Line Loop(ll2) = {-lr1,lr4,lr5,lr6};
	s2=news; Plane Surface(s2) = {ll2};
	Transfinite Surface{s2};
	Physical Surface("RB") = {s2};

	//---------Bottom----------------------
	lr7=newl; Line(lr7) = {p3,p31};
	lr8=newl; Line(lr8) = {p31,p23};

	Transfinite Line{lr7} = NbTrans;

	ll3=newll; Line Loop(ll3) = {-lr4,-l3,lr7,lr8};
	s3=news; Plane Surface(s3) = {ll3};
	Transfinite Surface{s3};
	Physical Surface("B") = {s3};

	//-----Bottom Left-----------------
	lr9=newl;  Line(lr9)  = {p3,p33};
	lr10=newl; Line(lr10) = {p33,p32};
	lr11=newl; Line(lr11) = {p32,p31};

	Transfinite Line{lr9,lr10,lr11}=NbTrans;

	ll4=newll; Line Loop(ll4) = {-lr7,lr9,lr10,lr11};
	s4=news; Plane Surface(s4) = {ll4};
	Transfinite Surface{s4};
	Physical Surface("LB") = {s4};

	// -------------Left---------------------------
	lr12=newl; Line(lr12) = {p4,p41};
	lr13=newl; Line(lr13) = {p41,p33};

	Transfinite Line{lr12}=NbTrans;

	ll5=newll; Line Loop(ll5) = {-lr9,-l2, lr12,lr13};
	s5=news; Plane Surface(s5) = {ll5};
	Transfinite Surface{s5};
	Physical Surface("L") = {s5};

	//-------Top Left-------------------
	lr14=newl; Line(lr14) = {p4,p43};
	lr15=newl; Line(lr15) = {p43,p42};
	lr16=newl; Line(lr16) = {p42,p41};

	Transfinite Line{lr14, lr15, lr16}=NbTrans;
	ll6=newll; Line Loop(ll6) = {lr14,lr15,lr16,-lr12};
	s6=news; Plane Surface(s6) = {ll6};
	Transfinite Surface{s6};
	Physical Surface("LT") = {s6};

	//------------Top-----------------------------
	lr17=newl; Line(lr17) = {p1,p13};
	lr18=newl; Line(lr18) = {p13,p43};

	Transfinite Line{lr17}=NbTrans;
	ll7=newll; Line Loop(ll7) = {lr17,lr18,-lr14,-l1};	
	s7=news; Plane Surface(s7) = {ll7};
	Transfinite Surface{s7};
	Physical Surface("T") = {s7};

	//------------Top Right-----------------------------
	lr19=newl; Line(lr19) = {p11,p12};
	lr20=newl; Line(lr20) = {p12,p13};

	Transfinite Line{lr19, lr20}=NbTrans;
	ll8=newll; Line Loop(ll8) = {lr19,lr20,-lr17,-lr3};	
	s8=news; Plane Surface(s8) = {ll8};
	Transfinite Surface{s8};
	Physical Surface("RT") = {s8};

	//
	//Transfinite Surface "*";
	//Physical Surface("PML") = {s1,s2,s3,s4,s5};
	// Bords du domaine PML absorbants
	Physical Line ("Left_PML")   = {lr10,lr13,lr16};
	Physical Line ("Right_PML")  = {lr2,lr6,lr19};
	Physical Line ("Bottom_PML") = {lr5,lr8,lr11};
	Physical Line ("Top_PML")    = {lr15,lr18,lr20};
EndIf
