
(* Dimitri Komatitsch, University of Pau, Dec 2007 *)

(* 'Mathematica' script to compute the spatial derivative of a Ricker in order
to calculate the analytical traction for a plane wave for Bielak's technique
on the absorbing edges *)

  rac3 = Sqrt[3];
  rac3sur2 = rac3 / 2;
  HALF = 1 / 2;

  rickertest[t_] := - (1 - 2 a t^2) Exp[-a t^2]

  Ux[t_,x_,z_] := rac3sur2 rickertest[t - x/2 + (9 - z) rac3sur2] + rac3sur2 rickertest[t - x/2 - (9 - z) rac3sur2] + rac3 rickertest[t - x/2]

  Uz[t_,x_,z_] := - HALF rickertest[t - x/2 + (9 - z) rac3sur2] + HALF rickertest[t - x/2 - (9 - z) rac3sur2]

(* check displacements Ux and Uz *)
(* Print[Simplify[Ux[t,x,z]]] *)

(* Print[Simplify[Uz[t,x,z]]] *)

(* output results to a file *)

" " >> result.txt

(* compute spatial derivatives of displacement *)

" " >>> result.txt
" dxUx = " >>> result.txt
" " >>> result.txt
 FortranForm[Simplify[D[Ux[t,x,z],x]]] >>> result.txt

" " >>> result.txt
" dzUx = " >>> result.txt
" " >>> result.txt
 FortranForm[Simplify[D[Ux[t,x,z],z]]] >>> result.txt

" " >>> result.txt
" dxUz = " >>> result.txt
" " >>> result.txt
 FortranForm[Simplify[D[Uz[t,x,z],x]]] >>> result.txt

" " >>> result.txt
" dzUz = " >>> result.txt
" " >>> result.txt
 FortranForm[Simplify[D[Uz[t,x,z],z]]] >>> result.txt

" " >>> result.txt

