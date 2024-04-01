(* ::Package:: *)

(* ::Input:: *)
(*(**)
(**)
(*\[Psi]o = fo[n] BesselJ[Abs[n],ko r];*)
(*\[Psi][s] = Ao[n]HankelH1[Abs[n],ko r];*)
(*\[Psi][0] = f[0,n] BesselJ[Abs[n],k[0] r];*)
(*\[Psi][1] = f[1,n] BesselJ[Abs[n],k[1] r] + A[1,n]HankelH1[Abs[n],k[1] r];*)
(**)
(*eqs = {\[Psi][1]  == \[Psi][0]/.r-> a[0],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][0],r]/\[Rho][0]/.r-> a[0], \[Psi][s] + \[Psi]o  == \[Psi][1]/.r-> a[1],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][s] + \[Psi]o ,r]/\[Rho]o/.r-> a[1]};*)
(**)
(*sol = Ao[n]/.Solve[eqs,{Ao[n],f[0,n],f[1,n],A[1,n]}] /.{\[Rho]o -> qo ko,\[Rho][1] -> q[1] k[1],\[Rho][0] -> q[0] k[0]} //Simplify;*)*)


(* ::Input:: *)
(*ClearAll[\[Psi]o,\[Psi],J,H]*)
(*\[Psi]o = fo[n] Subscript[J, n][ko r];*)
(*\[Psi][s] = Ao[n]Subscript[H, n][ko r];*)
(*\[Psi][0] = f[0,n] Subscript[J, n][k[0] r];*)
(*\[Psi][1] = f[1,n] Subscript[J, n][k[1] r] + A[1,n]Subscript[H, n][k[1] r];*)
(*eqs = {\[Psi][1]  == \[Psi][0]/.r-> a[0],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][0],r]/\[Rho][0]/.r-> a[0], \[Psi][s] + \[Psi]o  == \[Psi][1]/.r-> a[1],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][s] + \[Psi]o ,r]/\[Rho]o/.r-> a[1]};*)
(**)
(*subsol = Solve[eqs,{Ao[n],f[0,n],f[1,n],A[1,n]}] /.{\[Rho]o -> qo ko,\[Rho][1] -> q[1] k[1],\[Rho][0] -> q[0] k[0]} //Simplify;*)


(* ::Input:: *)
(*(*The scattered wave*)*)
(*Ao[n]/.subsol ;*)
(*tmp = Apart@FullSimplify[%[[1]]/.{Inactive[BesselJ]:> J,Inactive[HankelH1]:> H, a[n_]:> Subscript[a, n], q[n_]:> Subscript[q, n], k[n_]:> Subscript[k, n]}];*)


(* ::Input:: *)
(*F[f_,g_][x_,y_]:= f[x] g[y]-f[y] g[x]*)
(*Fd[f_,g_][x_,y_]:= f[x] g'[y]-f'[y] g[x]*)
(**)
(*Yd[x_,y_]=  Fd[Subscript[H, n],Subscript[J, n]] [x,y];*)
(*Y[x_,y_]=  F[Subscript[H, n],Subscript[J, n]] [x,y];*)
(*Ydd[x_,y_]=  F[Subscript[H, n]',Subscript[J, n]'] [x,y];*)
(*ClearAll[q]*)
(*dq = qo/Subscript[q, 1];*)
(*dq0 = Subscript[q, 0]/Subscript[q, 1];*)
(* dq/.{qo -> \[Rho]o co, Subscript[q, 1]-> Subscript[\[Rho], 1] Subscript[c, 1] , Subscript[q, 0]-> Subscript[\[Rho], 0] Subscript[c, 0]}*)
(*dq0/.{qo -> \[Rho]o co, Subscript[q, 1]-> Subscript[\[Rho], 1] Subscript[c, 1] , Subscript[q, 0]-> Subscript[\[Rho], 0] Subscript[c, 0]}*)
(**)
(*numer =Inactive[Yd][ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1],ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1]](-Inactive[Y][Subscript[a, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1],Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Derivative[1][Subscript[J, n]][Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 0]] + Subscript[q, 0]/Subscript[q, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[J, n][Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 0]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Inactive[Yd][Subscript[a, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1],Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1]])*)
(*denom = Derivative[1][Subscript[J, n]][Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 0]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)(Derivative[1][Subscript[H, n]][ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Inactive[Y][Subscript[a, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1],Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1]]+qo/Subscript[q, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[H, n][ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Inactive[Yd][Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1],Subscript[a, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1]])+Subscript[q, 0]/Subscript[q, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[J, n][Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 0]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)(qo/Subscript[q, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[H, n][ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Inactive[Ydd][Subscript[a, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1],Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1]] - Derivative[1][Subscript[H, n]][ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1]]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Inactive[Yd][Subscript[a, 1]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1],Subscript[a, 0]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[k, 1]])*)
(**)
(*tmp2 = Inactive[fo][n]\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)( numer -denom  Subscript[J, n][ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1]])/(denom Subscript[H, n][ko\!\(\**)
(*TagBox["*",*)
(*"InactiveToken",*)
(*BaseStyle->"Inactive",*)
(*SyntaxForm->"*"]\)Subscript[a, 1]])*)
(**)
(*tmp ==tmp2 //Activate//Simplify*)
(**)


(* ::Input:: *)
(*Activate[tmp2]/.{Subscript[q, 0] -> Subscript[q, 1],Subscript[k, 1] -> kin,Subscript[k, 0]-> kin} //FullSimplify*)


(* ::Input:: *)
(*(**)
(*eqs = {\[Psi][1]  == \[Psi][0]/.r-> a[0],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][0],r]/\[Rho][0]/.r-> a[0], \[Psi][s] + \[Psi]o  == \[Psi][1]/.r-> a[1], D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][s] + \[Psi]o ,r]/\[Rho]o/.r-> a[1]};*)
(*subsol = Solve[eqs,{Ao[n],f[0,n],f[1,n],A[1,n]}] /.{\[Rho]o -> qo ko,\[Rho][1] -> q[1] k[1],\[Rho][0] -> q[0] k[0]}//Simplify;*)
(**)*)
(**)
(*f[0,n]/.subsol [[1]]//Simplify;*)
(*f[1,n]/.subsol [[1]]//Simplify;*)
(*A[1,n]/.subsol [[1]]//Simplify;*)


(* ::Input:: *)
(*eqs/.subsol//.{ qo-> \[Rho]o /ko ,  q[1]->\[Rho][1]/k[1] ,q[0] -> \[Rho][0]/k[0]}//Simplify*)


(* ::Input:: *)
(*denom = dq0 Subscript[J, n][a[0] k[0]] (dq Subscript[J, n][ko a[1]] Inactive[Ydd][a[1] k[1],a[0] k[1]] - Derivative[1][Subscript[J, n]][ko a[1]]Inactive[Yd][a[1] k[1],a[0] k[1]] )+Derivative[1][Subscript[J, n]][a[0] k[0]](dq Subscript[J, n][ko a[1]] Inactive[Yd][a[0] k[1],a[1] k[1]]+ Derivative[1][Subscript[J, n]][ko a[1]]Inactive[Y][a[1] k[1],a[0] k[1]]) ;*)
(*force = Ao[n]  Yd[ko a[1],ko a[1]]/denom;*)
(**)
(*numer =- dq0  Yd[k[1] a[0],k[1] a[0]];*)
(*numer force == f[0,n]/.subsol2[[1]]/.{q[0] -> dq0 q[1],qo -> dq q[1] }//Activate//Simplify*)
(**)
(*numer =  dq0 Subscript[J, n][a[0] k[0]] Derivative[1][Subscript[H, n]][a[0] k[1]]-Subscript[H, n][a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]];*)
(*numer force == f[1,n]/.subsol2[[1]]/.{q[0] -> dq0 q[1],qo -> dq q[1] }//Activate//FullSimplify*)
(**)
(*numer =  Subscript[J, n][a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]]-dq0 Subscript[J, n][a[0] k[0]] Derivative[1][Subscript[J, n]][a[0] k[1]];*)
(*numer force == A[1,n]/.subsol2[[1]]/.{q[0] -> dq0 q[1],qo -> dq q[1] }//Activate//FullSimplify*)
(**)
(**)
(*A[1,n]/.subsol[[1]]//.{q[0] -> q[1],k[0]-> k[1],qo -> dq q[1] }//Activate//FullSimplify*)
(*f[1,n]/.subsol[[1]]//.{q[0] -> q[1],k[0]-> k[1],qo -> dq q[1] }//Activate//FullSimplify*)
(*f[0,n]/.subsol[[1]]//.{q[0] -> q[1],k[0]-> k[1],qo -> dq q[1] }//Activate//FullSimplify*)
(**)
(**)


(* ::Subsection:: *)
(*Capsule filled with void / air*)


(* ::Input:: *)
(*ClearAll[\[Psi]o,\[Psi],J,H]*)
(*\[Psi]o = fo[n] Subscript[J, n][ko r];*)
(*\[Psi][s] = Ao[n]Subscript[H, n][ko r];*)
(*\[Psi][1] = f[1,n] Subscript[J, n][k[1] r] + A[1,n]Subscript[H, n][k[1] r];*)
(*eqs = {\[Psi][1]  /.r-> a[0], \[Psi][s] + \[Psi]o  - \[Psi][1]/.r-> a[1],D[\[Psi][1],r]/\[Rho][1]  - D[\[Psi][s] + \[Psi]o ,r]/\[Rho]o/.r-> a[1]};*)
(**)
(*subsol = Solve[Thread[eqs==0],{Ao[n],f[1,n],A[1,n]}]  //Simplify;*)


(* ::Input:: *)
(*subN = {a[1]->2.0,a[0]->1.5,\[Rho]o->0.00001,\[Rho][1]->2.0,ko-> 0.0001,k[1] -> 2.0,Subscript[J, n][ko a[1]]->1.0,fo[n]->1, Subscript[H, n_][x_] ->HankelH1[n,x], Subscript[J, n_][x_] ->BesselJ[n,x], Derivative[1][Subscript[H, n_]][x_] ->D[HankelH1[n,x],x],Derivative[1][Subscript[J, n_]][x_] ->D[BesselJ[n,x],x]};*)
(*subsol//.subN//Flatten;*)


(* ::Input:: *)
(*%/.n->#&/@Range[1,20,3]*)
(**)


(* ::Input:: *)
(*(*Capsule with given forcing on boundary*)*)
(*ClearAll[\[Psi]o,\[Psi],J,H]*)
(*\[Psi]o = fo[n] Subscript[J, n][ko r];*)
(*\[Psi][1] = f[1,n] Subscript[J, n][k[1] r] + A[1,n]Subscript[H, n][k[1] r];*)
(*eqs = {\[Psi][1]  /.r-> a[0],  \[Psi]o  - \[Psi][1]/.r-> a[1]};*)
(*subsol = Solve[Thread[eqs==0],{f[1,n],A[1,n]}] /.{\[Rho]o -> qo ko,\[Rho][1] -> q[1] k[1],\[Rho][0] -> q[0] k[0]} //Simplify*)
(**)
(*vars = {f[1,n],A[1,n]};*)
(*M = Coefficient[eqs,#]&/@vars;*)
(*M = Transpose@M;*)
(*b =  eqs - M . vars //Simplify;*)
(*M//MatrixForm*)


(* ::Input:: *)
(*subN = {a[1]->2.0,a[0]->1.5,k[1] -> 2.0,Subscript[J, n][ko a[1]]->1.0,fo[n]->1, Subscript[H, n_][x_] ->HankelH1[n,x], Subscript[J, n_][x_] ->BesselJ[n,x]};*)
(*ns = Range[1,30,3];*)
(*subsol//.subN//Flatten;*)
(*subNsol=Flatten[%/.n->#&/@ns];*)
(**)
(*eqs//.subN/.n->#&/@ns;*)
(*%/.subNsol*)


(* ::Input:: *)
(**)
(*{Abs@f[1,n]Abs@BesselJ[n,a[1] k[1]],Abs@A[1,n]Abs@HankelH1[n,a[1] k[1]]}//.subN//Flatten;*)
(*%/.n->#&/@ns;*)
(*%/.Flatten@subNsol*)
(**)


(* ::Input:: *)
(*(*Is the matrix system ill-posed?*)*)
(*NM = M //.subN;*)
(*Nb = b //.subN;*)
(*ns = Range[1,70,6];*)
(*{Abs@Det[NM/.n->#],SingularValueList[NM/.n->#]}&/@ns*)
(**)
(*(*If we numerically solve the system:*)*)
(*sols =Flatten[ Thread[(vars/.subN/.n->#)->LinearSolve[NM/.n->#,Nb/.n->#]]&/@ns];*)
(**)
(**)


(* ::Input:: *)
(*eqs//.subN/.n->#&/@ns;*)
(*%/.sols*)
(**)


(* ::Input:: *)
(*subsol*)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


Test a capsule filled with air





(* ::Input:: *)
(*(* For one particle *)*)
(*(*Clear[dq0,dq]*)
(*denom = Y[a[1] k[1],a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]]-dq0 Subscript[J, n][a[0] k[0]] Yd[a[1] k[1],a[0] k[1]];*)
(*force = (Ao[n] Subscript[H, n][ko a[1]]+fo[n] Subscript[J, n][ko a[1]])/denom;*)
(**)
(*numer = - dq0  Yd[a[0] k[1],a[0] k[1]]*)
(*force numer  == f[0,n]/.subsol [[1]]/.{dq0 \[Rule] q[0]/q[1]}//Simplify*)
(*numer = dq0 Subscript[J, n][a[0] k[0]] Derivative[1][Subscript[H, n]][a[0] k[1]]- Subscript[H, n][a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]];*)
(*force numer  ==f[1,n]/.subsol [[1]]/.{dq0 \[Rule] q[0]/q[1]}//Simplify*)
(**)
(*numer =  Subscript[J, n][a[0] k[1]] Derivative[1][Subscript[J, n]][a[0] k[0]]-dq0 Subscript[J, n][a[0] k[0]] Derivative[1][Subscript[J, n]][a[0] k[1]]*)
(*force  numer== A[1,n]/.subsol [[1]]/.{dq0 \[Rule] q[0]/q[1]}//Simplify*)
(**)*)


(* ::Input:: *)
(*D[Y[x,y],x,y]*)
(**)
(*Yd [y,x]+Yd [x,y]//Simplify*)
(**)
(*capsule = tmp2;*)


Test subtracting the frequency response of the capsule . 
I'll assume the background is water, that there is Sunflower oil
in the tube, and the tube itself is silica;




subN = {Subscript[a, 1] -> 1, Subscript[a, 0] -> 3/4., qo-> 1, Subscript[q, 1]-> 6, Subscript[q, 0]-> 0.9, Subscript[k, 1]-> ko/3., Subscript[k, 0] -> 1.1 ko, fo[_] -> 1};
subHJ = { 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{
RowBox[{"Derivative", "[", "1", "]"}], "[", 
RowBox[{"Subscript", "[", 
RowBox[{"H", ",", "n_"}], "]"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) -> D[HankelH1[n,\[FormalK]],\[FormalK]]/.\[FormalK]->k, 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{
RowBox[{"Derivative", "[", "1", "]"}], "[", 
RowBox[{"Subscript", "[", 
RowBox[{"J", ",", "n_"}], "]"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) ->  D[BesselJ[n,\[FormalK]],\[FormalK]]/.\[FormalK]->k, 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{"Subscript", "[", 
RowBox[{"H", ",", "n_"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) :> HankelH1[n,k], 
\!\(\*
TagBox[
StyleBox[
RowBox[{
RowBox[{"Subscript", "[", 
RowBox[{"J", ",", "n_"}], "]"}], "[", "k_", "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\) :> BesselJ[n,k]
};

tube = capsule/.{Subscript[k, 0]-> ko, Subscript[q, 0]-> qo} //FullSimplify;
core = capsule/.{Subscript[k, 1]-> ko, Subscript[q, 1]-> qo} //FullSimplify;

tubeSubtract = capsule - tube  /.subHJ/.subN //Simplify;
Ncore = core/.subHJ/.subN //Simplify;


subN = {Subscript[a, 1] -> 1, Subscript[a, 0] -> 2., qo-> 4., 
    Subscript[q, 1]-> 1, Subscript[q, 0]-> 9., ko->0.5,
   Subscript[k, 1]-> 1., Subscript[k, 0] -> 1/3., fo[_] -> 1
};

Activate@capsule/.subHJ/.subN/.{n->1}




Activate@Inactive[Y][1\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`,2.`\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`]/.subHJ/.n->1
Activate@Ydd[1\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`,2.`\!\(\*
TagBox["*",
"InactiveToken",
BaseStyle->"Inactive",
SyntaxForm->"*"]\)1.`]/.subHJ/.n->1


Test subtracting the frequency response of the capsule . 
I'll assume the background is water, that there is Sunflower oil
in the tube, and the tube itself is silica;

Plot[{Abs@Ncore/.ko-> k/.n-> 0,  Abs@tubeSubtract/.ko-> k/.n-> 0},{k,0,1.}, 
PlotLegends->{"Just core", "Capsule - Tube"}, AxesLabel->{"k x tube radius", "response"} ]
Legended[Graphics[{{{}, {}, Annotation[{RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["


(* Rationale:  

Let u[c,t] be the reponse of core + tube (the capsule), u[0,t] be the reponse of just the tube,
and u[c,0] be the reponse of just the core, so that u[0,0] =0 
(u[0,0] the response with the tube and core being the same as the background material.

Then 
(1) u[c,t] - u[0,t] = \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, t]\)\)c + O[c^2] = \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, 0]\)\)c + \!\(
\*SubscriptBox[\(\[PartialD]\), \(ct\)]\(u[0, 0]\)\)c t + O[c t] + O[c^2]
(2) u[c,0] = u[0,0] + \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, 0]\)\)c + O[c^2] = \!\(
\*SubscriptBox[\(\[PartialD]\), \(c\)]\(u[0, 0]\)\)c + O[c^2]
using (1) and (2) implies that 
u[c,t] - u[0,t] = O[c t] + O[c^2]
So the two are similar if the material in the tube is similar to 
the background material.  
    















(* ::Section:: *)
(*Addition Theorems*)


(* ::Input:: *)
(*ClearAll[Graf,R,\[CapitalTheta]]*)
(*R[l_] := Norm[{x,y}-{x[l],y[l]}];*)
(*\[CapitalTheta][l_] := ArcTan@@({x,y}-{x[l],y[l]});*)
(*R[l_,j_] := Norm[{x[j],y[j]}-{x[l],y[l]}];*)
(*\[CapitalTheta][l_,j_] := ArcTan@@({x[j],y[j]}-{x[l],y[l]});*)
(**)
(*Graf[H_,M_]:= H[n,R[l]]E^(I n \[CapitalTheta][l])  - Sum[H[n-m,R[l,j]]E^(I (n-m) \[CapitalTheta][l,j]) BesselJ[m, R[j]]E^(I m \[CapitalTheta][j]),{m,-M,M}];*)
(*error =Abs[ Graf[HankelH1,50]/HankelH1[n,R[l]] ]/.{x[l]-> 0,y[l]-> 0, x[j]-> 2,y[j]-> 0, y->0};*)
(*Plot[error/.{n->3},{x,0.1,1},PlotRange-> All]*)
(**)
(*error =Abs[ Graf[HankelH1,3]/HankelH1[n,R[l]] ]/.{x[l]-> 0,y[l]-> 0, x[j]-> 2,y[j]-> 0, y->0};*)
(*Plot[error/.{n->3},{x,1.,2},PlotRange-> All]*)
(**)


(* ::Input:: *)
(*error =Abs[ Graf[BesselJ,20]/BesselJ[n,R[l]] ]/.{x[l]-> 0,y[l]-> 0, x[j]-> 2,y[j]-> 0, y->0};*)
(*Plot[error/.{n->3},{x,0.3,2},PlotRange-> All]*)
(**)
(**)
