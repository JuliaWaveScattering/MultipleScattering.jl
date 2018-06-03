(* ::Package:: *)

(* ::Input:: *)
(*$Assumptions ={ao>0,a[_]>0,ko>0,k[_]>0,c[_]>0,co>0, n \[Element] Integers};*)
(**)
(*\[Psi]o = fo[n] BesselJ[Abs[n],ko r];*)
(*\[Psi][s] = Ao[n]HankelH1[Abs[n],ko r];*)
(*\[Psi][0] = f[0,n] BesselJ[Abs[n],k[0] r];*)
(*\[Psi][1] = f[1,n] BesselJ[Abs[n],k[1] r] + A[1,n]HankelH1[Abs[n],k[1] r];*)
(**)
(**)
(*eqs = {\[Psi][1]  == \[Psi][0]/.r-> a[0],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][0],r]/\[Rho][0]/.r-> a[0], \[Psi][s] + \[Psi]o  == \[Psi][1]/.r-> a[1],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][s] + \[Psi]o ,r]/\[Rho]o/.r-> a[1]};*)
(**)
(*sol = Ao[n]/.Solve[eqs,{Ao[n],f[0,n],f[1,n],A[1,n]}] /.{\[Rho]o -> qo ko,\[Rho][1] -> q[1] k[1],\[Rho][0] -> q[0] k[0]} //Simplify;*)


(* ::Input:: *)
(*ClearAll[\[Psi]o,\[Psi],J,H]*)
(*\[Psi]o = fo[n] Subscript[J, n][ko r];*)
(*\[Psi][s] = Ao[n]Subscript[H, n][ko r];*)
(*\[Psi][0] = f[0,n] Subscript[J, n][k[0] r];*)
(*\[Psi][1] = f[1,n] Subscript[J, n][k[1] r] + A[1,n]Subscript[H, n][k[1] r];*)
(*eqs = {\[Psi][1]  == \[Psi][0]/.r-> a[0],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][0],r]/\[Rho][0]/.r-> a[0], \[Psi][s] + \[Psi]o  == \[Psi][1]/.r-> a[1],D[\[Psi][1],r]/\[Rho][1]  == D[\[Psi][s] + \[Psi]o ,r]/\[Rho]o/.r-> a[1]};*)
(*Ao[n]/.Solve[eqs,{Ao[n],f[0,n],f[1,n],A[1,n]}] /.{\[Rho]o -> qo ko,\[Rho][1] -> q[1] k[1],\[Rho][0] -> q[0] k[0]} ;*)
(**)
(*tmp = Apart@FullSimplify[%[[1]]/.{Inactive[BesselJ]:> J,Inactive[HankelH1]:> H, a[n_]:> Subscript[a, n], q[n_]:> Subscript[q, n], k[n_]:> Subscript[k, n]}]*)


(* ::Input:: *)
(*tmp == fo[n]/Subscript[H, n][ko Subscript[a, 1]] (-Subscript[J, n][ko Subscript[a, 1]]+Subscript[q, 1]((Subscript[J, n][ko Subscript[a, 1]] Derivative[1][Subscript[H, n]][ko Subscript[a, 1]]- Derivative[1][Subscript[J, n]][ko Subscript[a, 1]]Subscript[H, n][ko Subscript[a, 1]]) (Subscript[q, 1] (Subscript[H, n][Subscript[a, 1] Subscript[k, 1]] Subscript[J, n][Subscript[a, 0] Subscript[k, 1]]-Subscript[H, n][Subscript[a, 0] Subscript[k, 1]] Subscript[J, n][Subscript[a, 1] Subscript[k, 1]]) Derivative[1][Subscript[J, n]][Subscript[a, 0] Subscript[k, 0]]+Subscript[q, 0] Subscript[J, n][Subscript[a, 0] Subscript[k, 0]] (Subscript[J, n][Subscript[a, 1] Subscript[k, 1]] Derivative[1][Subscript[H, n]][Subscript[a, 0] Subscript[k, 1]]-Subscript[H, n][Subscript[a, 1] Subscript[k, 1]] Derivative[1][Subscript[J, n]][Subscript[a, 0] Subscript[k, 1]])) )/(Subscript[q, 1] Derivative[1][Subscript[J, n]][Subscript[a, 0] Subscript[k, 0]] (Subscript[q, 1] Derivative[1][Subscript[H, n]][ko Subscript[a, 1]](Subscript[H, n][Subscript[a, 1] Subscript[k, 1]] Subscript[J, n][Subscript[a, 0] Subscript[k, 1]]-Subscript[H, n][Subscript[a, 0] Subscript[k, 1]] Subscript[J, n][Subscript[a, 1] Subscript[k, 1]]) +qo Subscript[H, n][ko Subscript[a, 1]] (Subscript[H, n][Subscript[a, 0] Subscript[k, 1]] Derivative[1][Subscript[J, n]][Subscript[a, 1] Subscript[k, 1]]-Subscript[J, n][Subscript[a, 0] Subscript[k, 1]] Derivative[1][Subscript[H, n]][Subscript[a, 1] Subscript[k, 1]]))+Subscript[q, 0] Subscript[J, n][Subscript[a, 0] Subscript[k, 0]] (Subscript[q, 1] Derivative[1][Subscript[H, n]][ko Subscript[a, 1]] (Subscript[J, n][Subscript[a, 1] Subscript[k, 1]] Derivative[1][Subscript[H, n]][Subscript[a, 0] Subscript[k, 1]]-Subscript[H, n][Subscript[a, 1] Subscript[k, 1]] Derivative[1][Subscript[J, n]][Subscript[a, 0] Subscript[k, 1]])+qo Subscript[H, n][ko Subscript[a, 1]] (Derivative[1][Subscript[H, n]][Subscript[a, 1] Subscript[k, 1]] Derivative[1][Subscript[J, n]][Subscript[a, 0] Subscript[k, 1]]-Derivative[1][Subscript[H, n]][Subscript[a, 0] Subscript[k, 1]] Derivative[1][Subscript[J, n]][Subscript[a, 1] Subscript[k, 1]]))))//Simplify*)


(* ::Input:: *)
(*F[f_,g_][x_,y_]:= f[x] g[y]-f[y] g[x]*)
(*Fd[f_,g_][x_,y_]:= f[x] g'[y]-f'[y] g[x]*)
(**)
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
(*tmp == -fo[n] Subscript[J, n][ko Subscript[a, 1]]/Subscript[H, n][ko Subscript[a, 1]]-fo[n] Yd [ko Subscript[a, 1],ko Subscript[a, 1]] /Subscript[H, n][ko Subscript[a, 1]] (  Y [Subscript[a, 1] Subscript[k, 1],Subscript[a, 0] Subscript[k, 1]] Derivative[1][Subscript[J, n]][Subscript[a, 0] Subscript[k, 0]]-dq0 Subscript[J, n][Subscript[a, 0] Subscript[k, 0]] Yd [Subscript[a, 1] Subscript[k, 1],Subscript[a, 0] Subscript[k, 1]] )/( Derivative[1][Subscript[J, n]][Subscript[a, 0] Subscript[k, 0]] (dq Subscript[H, n][ko Subscript[a, 1]] Yd [Subscript[a, 0] Subscript[k, 1],Subscript[a, 1] Subscript[k, 1]]+ Derivative[1][Subscript[H, n]][ko Subscript[a, 1]]Y[Subscript[a, 1] Subscript[k, 1],Subscript[a, 0] Subscript[k, 1]])+dq0 Subscript[J, n][Subscript[a, 0] Subscript[k, 0]] (dq Subscript[H, n][ko Subscript[a, 1]]Ydd [Subscript[a, 1] Subscript[k, 1],Subscript[a, 0] Subscript[k, 1]]- Derivative[1][Subscript[H, n]][ko Subscript[a, 1]] Yd[Subscript[a, 1] Subscript[k, 1],Subscript[a, 0] Subscript[k, 1]])) //Simplify*)
(**)


(* ::Input:: *)
(*D[Y[x,y],x,y]*)
(**)
(*Yd [y,x]+Yd [x,y]//Simplify*)
(**)
(*capsule = tmp;*)


Test subtracting the frequency response of the capsule. 
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



Test subtracting the frequency response of the capsule. 
I'll assume the background is water, that there is Sunflower oil
in the tube, and the tube itself is silica;

Plot[{Abs@Ncore/.ko-> k/.n-> 0,  Abs@tubeSubtract/.ko-> k/.n-> 0},{k,0,1.}, 
PlotLegends->{"Just core", "Capsule - Tube"}, AxesLabel->{"k x tube radius", "response"} ]
Legended[Graphics[{{{}, {}, Annotation[{RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["
1:eJwVy3c81QsDx3GrEmWlRCHXqJAkpPX7iqRQNgkHGaESlYZuSRHKKE5GZlZG
OpzlGkVWKTsjW4e6CD8NNMRznz8+r/dfH7mT5yzcebi4uFL+6/8auY+1VY1b
E5Vbj6mqxITu49OX8h2S1ITTiL1dQKc+sTMtIZkheRAG0zv+yHe6EqMbFYIu
Kx2DXbS1Qa7pJcK4qZueImkFhQdft//suE0oREhnL8naQaD5q9k/C1HEgpFr
vIsSBSJu4v1hpnHERHBHrUucC3JfCbTyUFIJW7f55zGSbtAeeCr2sSOTqNWX
KqlJ9gDDSvfmOck8Ql1+f9F3WS/oej2pTV94RiTzOOcpZp6GMKswYF8lneDn
3MqwUfJB/PM88b9N2UTZcHtTjZUvqGZX9wlvLyN8zdyMHOL84FD0e3HW8QUx
sD1UP1LyIh4seRgMddQQ1LT1VYr2/hgXvalj5lpPGAnn73uRfAm3rjavPS3Z
QLCnG7VI2asI9FlJS11oJs5QKPRQlwBcXdb8vE+8jfirmVTblHkNabGdW/dX
thNRhWJbzJVu4EuTlGKgaRdhIJ2VNX4qEHSJ1OdOWu+JG/WPRF5a3YRFm6fP
2u29hKa2aqRdXBDOM+6oclEGifS1pnfuSgbDLWFV53jHCHGtKb7WpiwYyody
G9QkPhI2IcM88vYhKL1h5ERx/USsmvULLE++gw0eg1euS44TVzqoVyZlw9Ej
9i+rfGGasIwYYP9TFQ7llHBuFcMZQu2g0mywy12w35sFT4l/IUYZJX7SmfdQ
sP5hlkXlV8I8ttf7mFIUjA49vpZiOkeoGsvnSb2KwoY3FrS1uvMEP++Zfz+d
ioZ5ZJPxba0fhFPpYnGl1X20Gj6uI7b/ItZuva92PO4BdnK6JrQpi0QQP2Nz
uCQVapkOpes6eaFs+CCM7kuF1E6yNOoWH9pDzo33vaIiIVagZkZjGRR4VQvU
Lj1ETv/Bb+eTl6P+T5ZaR3scinxqJs7FrMSqb3FashGJ4PicjFdaJgLWDv+4
wyOJiMi5LLPzgggovpbzfnsewfUjnc7kiIA2JVxaO/YI3b/b+STqRWExFrbP
2yAZhWs0nI6nrkHCQIA+azEVc6T2qWaqBA5stMsatE7DJiKJR1VqPSZO7FrG
X5iG1yFDjMX09ZCZvOpauC0drWs3n+BlSOKcl7YiV386xFy8a8PHNkDMjZab
rZMBAdv8l3Umm3Dc7jFt5ksWVrZ8FmKEK4LssxLIkMqGoDy/eMmAIu5Q+N0t
9bMhO6Nb46ihBIbrOSk2NRtPfrA3xQ4oQchnf8i1XTmQ9XLrT9m7BXW3ek4s
v/4ElwRZATvXqUDjqejyjfz50BQWUC3KUEeDWp1zk3o+bI/JzfV3qcO5+Er5
Dbt8JGSFStgL7kAke8j3Q34+GKOXnI39d2CsqrAv52gBtM7or1ow0UBq55Hi
HbFPQTTldKYIakJwKcjBUIaGuaJin5G6XUjwsqWoHqJB06/FgZtXBwodqs6i
PjRIPDUJOKurAyKv27XvOQ2mNj766hU6KM2LY5WoFcFo8xUh55LdmLKcrDBu
LkJuVX2oZ8VeWOfGv72wmo4gTmBAFJcukh1vmXjI0UGv2+XBpaALjtjZpuNa
dMgXe3olGuri3HW9lv2OdMT0yrF1o3URbjHVvqKQji159p5umw7g+YJeb5IJ
A/zLY5S5jPSgaD49VnOPicNF8jNhZQchFWi459RjJjx4hLOMPh2ESGH6PYES
JvoFo2z2iBngF7+lmjmHCdc2gZJMbwO0VpVcGNRhYa/TIbdU6UMIUA9a/PGR
hY7f35x3RhiiWVhcXE23BAlfH9j33TPC5ebdSPheCr36pcWou2YQUmN9bOAv
Q+/DkSybLDNkRapH/N5YBn/zMlONF2ZoM9ncQzEog5gMH1PxqxlU3opfUIwr
g8jKE6yFE+YYekXm0HeV45WjTORpdQsYvswWagqoQPXKIUmTCUtIMEUHubkr
wXyT4Sf+yBZsjn5jy75qPJx24Gywdka1d8oGu8la7DWp2hx6wx2zSvuFC8++
wu2hqtmGL97g+VkgskHxDRqr56S76nyhveTMtnzWiMShkQHvbRfRGKzZXK/c
gqSjB/fItl2G7P3PP/OmW6E35x4vJ/U3Ps8PDXaqt2PG4+NdQa2byJoofyMU
9Q65mrUFH9i3ICT9hBwd7UDf6TDjguMhiPBieaQrdOF6Ke+xkm+hmNMtGPUM
7EbXhZrzP2fu4qJGLXNdzXtMLr/Rs5oTiaFntffa5HqhpHjZYyrlPn57dyRt
DehDlFSK/bBtDFQOrFj0KO/HaKL/3bDlVBhmyJ0qlhhEVI6a+9aph7gs7fTX
B4chaEm87mfUxUN10lqw++gwTEVlK5jZiVgjPMNf+2sYJo6N8y1VSeCrTKLJ
J3/AgRXzL16/TsFpzfKelcYcVNDdypRWpKPgyMn0F3McqHfbMP23PcbsZymR
T3EjIMoS1+UZZoBlachLPTgKdldPdldQJgIPh6QxP43CWJ56LfRZFuQGjN5n
hnxELqNsPW9/NvbsCqZKanzCvVE+hVPKT+Ag0/KT+90n0OgxH2Y8czH7JpLj
G/Qvoi3PExINeUjgo1y0UBhDTGxH2ur1BVgzyMrgqR3D4Xzzwkn3p1h2Pchq
t984an/YiI10FoIrj0LjiE1gpviGzyUNGs4/okWKVU/gisl2vYrjRXiuEsYv
6vMZ5cvztPdJFEPU6LzqLaFJbJNOpSy9KUb/ttZ3Vysm8VKH6fgrmo6iVCiM
UqaQfEiDk32EAUeatqj04hTOyn3uV+Zm4mbI46D+/GlIGphG1TcwUTgmnSRk
QiI5ReqV5h0WwlJCnAePkSjZ05hyM5wFN4spxWfmJPyrv4c1RrCwsaKi6Kgt
iYxW32S3WBbuRtvXR5wksSvqJvthOgue2olfBK6SSB0/kf+7jIW/gsUPr8gh
kd/dqtJFsrCo8/fq7lwSb/Vc/BW+sdAzNdKeU0CiX+xR2/k5Fu7bMhwPFZNo
4n77UvgPC0vKFhdCKkjYT1MDjQTY6G+LTuV9R8Llt4xbjQIb/9yZd33XSaKX
/2yb2BY2qHudtma+J6E4es3spAobxtlqTL3B/36JjNvcGmyUXWlqCBon8WUm
oEIXbMRt04o2nyThMLWbfV+PDT9OspUcScK51+DlsAEbW0zODFV9J3GBKiAY
ZMIGH3dH1v15EuaeOUfaTNkYZu31dv5F4vAe60Q5SzYqvDO3q/8hcSKS/5ef
DRvxsoKzS0skPEaLT1fbsfE/QTmwgQ==
"]]}, "Charting`Private`Tag$2929#1"], Annotation[{RGBColor[0.880722, 0.611041, 0.142051], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["
1: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
ErvoPn+R1jykqS+enJv77bMaTRY58PA/1XezIA==
"]]}, "Charting`Private`Tag$2929#2"]}, {}, {}}, AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {"k x tube radius", "response"}, AxesOrigin -> {0, 0}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[Part[#, 1]], (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[Part[#, 1]], (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[Part[#, 2]]}& )}}, PlotRange -> {{0, 1.}, {0., 0.10664428602978178`}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Placed[\!\(TraditionalForm\`FormBox[TemplateBox[{"\<\"Just core\"\>", "\<\"Capsule - Tube\"\>"}, "\<LineLegend\>", DisplayFunction -> \((FormBox[StyleBox[StyleBox[PaneBox[TagBox[GridBox[{{TagBox[GridBox[{{GraphicsBox[{{Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[0]]], PointSize[0.5], Opacity[1. ], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6]], {LineBox[{{0, 10}, {20, 10}}]}}, {Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[0]]], PointSize[0.5], Opacity[1. ], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6]], {}}}, AspectRatio -> Full, ImageSize -> {20, 10}, PlotRangePadding -> None, ImagePadding -> Automatic, BaselinePosition -> \((Scaled[0.1] -> Baseline)\)], #}, {GraphicsBox[{{Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[0]]], PointSize[0.5], Opacity[1. ], RGBColor[0.880722, 0.611041, 0.142051], AbsoluteThickness[1.6]], {LineBox[{{0, 10}, {20, 10}}]}}, {Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[0]]], PointSize[0.5], Opacity[1. ], RGBColor[0.880722, 0.611041, 0.142051], AbsoluteThickness[1.6]], {}}}, AspectRatio -> Full, ImageSize -> {20, 10}, PlotRangePadding -> None, ImagePadding -> Automatic, BaselinePosition -> \((Scaled[0.1] -> Baseline)\)], #2}}, GridBoxAlignment -> {"\<Columns\>" -> {Center, Left}, "\<Rows\>" -> {{Baseline}}}, AutoDelete -> False, GridBoxDividers -> {"\<Columns\>" -> {{False}}, "\<Rows\>" -> {{False}}}, GridBoxItemSize -> {"\<Columns\>" -> {{All}}, "\<Rows\>" -> {{All}}}, GridBoxSpacings -> {"\<Columns\>" -> {{0.5}}, "\<Rows\>" -> {{0.8}}}], "\<Grid\>"]}}, GridBoxAlignment -> {"\<Columns\>" -> {{Left}}, "\<Rows\>" -> {{Top}}}, AutoDelete -> False, GridBoxItemSize -> {"\<Columns\>" -> {{Automatic}}, "\<Rows\>" -> {{Automatic}}}, GridBoxSpacings -> {"\<Columns\>" -> {{1}}, "\<Rows\>" -> {{0}}}], "\<Grid\>"], Alignment -> Left, AppearanceElements -> None, ImageMargins -> {{5, 5}, {5, 5}}, ImageSizeAction -> "\<ResizeToFit\>"], LineIndent -> 0, StripOnInput -> False], {FontFamily -> "\<Arial\>"}, Background -> Automatic, StripOnInput -> False], TraditionalForm] &)\), Editable -> True, InterpretationFunction :> \((RowBox[{"\<LineLegend\>", "\<[\>", RowBox[{RowBox[{"\<{\>", RowBox[{RowBox[{"\<Directive\>", "\<[\>", RowBox[{RowBox[{"\<Opacity\>", "\<[\>", "\<1.`\>", "\<]\>"}], "\<,\>", InterpretationBox[ButtonBox[TooltipBox[GraphicsBox[{{GrayLevel[0], RectangleBox[{0, 0}]}, {GrayLevel[0], RectangleBox[{1, \(-1\)}]}, {RGBColor[0.368417, 0.506779, 0.709798], RectangleBox[{0, \(-1\)}, {2, 1}]}}, AspectRatio -> 1, Frame -> True, FrameStyle -> RGBColor[0.24561133333333335`, 0.3378526666666667, 0.4731986666666667], FrameTicks -> None, PlotRangePadding -> None, ImageSize -> Dynamic[{Automatic, 1.35  CurrentValue["\<FontCapHeight\>"]/AbsoluteCurrentValue[Magnification]}]], "\<RGBColor[0.368417, 0.506779, 0.709798]\>"], Appearance -> None, BaseStyle -> {}, BaselinePosition -> Baseline, DefaultBaseStyle -> {}, ButtonFunction :> With[{Typeset`box$ = EvaluationBox[]}, If[Not[AbsoluteCurrentValue["\<Deployed\>"]], SelectionMove[Typeset`box$, All, Expression]; FrontEnd`Private`$ColorSelectorInitialAlpha = 1; FrontEnd`Private`$ColorSelectorInitialColor = RGBColor[0.368417, 0.506779, 0.709798]; FrontEnd`Private`$ColorSelectorUseMakeBoxes = True; MathLink`CallFrontEnd[FrontEnd`AttachCell[Typeset`box$, FrontEndResource["\<RGBColorValueSelector\>"], {0, {Left, Bottom}}, {Left, Top}, "\<ClosingActions\>" -> {"\<SelectionDeparture\>", "\<ParentChanged\>", "\<EvaluatorQuit\>"}]]]], BaseStyle -> Inherited, Evaluator -> Automatic, Method -> "\<Preemptive\>"], RGBColor[0.368417, 0.506779, 0.709798], Editable -> False, Selectable -> False], "\<,\>", RowBox[{"\<AbsoluteThickness\>", "\<[\>", "\<1.6`\>", "\<]\>"}]}], "\<]\>"}], "\<,\>", RowBox[{"\<Directive\>", "\<[\>", RowBox[{RowBox[{"\<Opacity\>", "\<[\>", "\<1.`\>", "\<]\>"}], "\<,\>", InterpretationBox[ButtonBox[TooltipBox[GraphicsBox[{{GrayLevel[0], RectangleBox[{0, 0}]}, {GrayLevel[0], RectangleBox[{1, \(-1\)}]}, {RGBColor[0.880722, 0.611041, 0.142051], RectangleBox[{0, \(-1\)}, {2, 1}]}}, AspectRatio -> 1, Frame -> True, FrameStyle -> RGBColor[0.587148, 0.40736066666666665`, 0.09470066666666668], FrameTicks -> None, PlotRangePadding -> None, ImageSize -> Dynamic[{Automatic, 1.35  CurrentValue["\<FontCapHeight\>"]/AbsoluteCurrentValue[Magnification]}]], "\<RGBColor[0.880722, 0.611041, 0.142051]\>"], Appearance -> None, BaseStyle -> {}, BaselinePosition -> Baseline, DefaultBaseStyle -> {}, ButtonFunction :> With[{Typeset`box$ = EvaluationBox[]}, If[Not[AbsoluteCurrentValue["\<Deployed\>"]], SelectionMove[Typeset`box$, All, Expression]; FrontEnd`Private`$ColorSelectorInitialAlpha = 1; FrontEnd`Private`$ColorSelectorInitialColor = RGBColor[0.880722, 0.611041, 0.142051]; FrontEnd`Private`$ColorSelectorUseMakeBoxes = True; MathLink`CallFrontEnd[FrontEnd`AttachCell[Typeset`box$, FrontEndResource["\<RGBColorValueSelector\>"], {0, {Left, Bottom}}, {Left, Top}, "\<ClosingActions\>" -> {"\<SelectionDeparture\>", "\<ParentChanged\>", "\<EvaluatorQuit\>"}]]]], BaseStyle -> Inherited, Evaluator -> Automatic, Method -> "\<Preemptive\>"], RGBColor[0.880722, 0.611041, 0.142051], Editable -> False, Selectable -> False], "\<,\>", RowBox[{"\<AbsoluteThickness\>", "\<[\>", "\<1.6`\>", "\<]\>"}]}], "\<]\>"}]}], "\<}\>"}], "\<,\>", RowBox[{"\<{\>", RowBox[{#, "\<,\>", #2}], "\<}\>"}], "\<,\>", RowBox[{"\<LegendMarkers\>", "\<->\>", "\<None\>"}], "\<,\>", RowBox[{"\<LabelStyle\>", "\<->\>", RowBox[{"\<{\>", "\<}\>"}]}], "\<,\>", RowBox[{"\<LegendLayout\>", "\<->\>", "\<\"Column\"\>"}]}], "\<]\>"}] &)\)], TraditionalForm]\),After]]

Plot[ Abs[tubeSubtract/Ncore - 1]/.ko-> k/.n-> 0,{k,0,1.}, 
PlotLegends->{"Relative error of (Capsule - Tube)"}, AxesLabel->{"k x tube radius", "response"}]

Legended[Graphics[{{{}, {}, Annotation[{RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["
1:eJwVx3k41IsCh3HSolSWkiX7dhOp7IX5Kmu2m6VQsjRGhSyRKMXIWMraESoi
S4U0E8YSUZlEh2QoMbYsE8JPyXFaDvfcP97n87zyJwMdKKt4eHji/+3/WlEm
u55PHSGZpycHGs9+wGoTyaBhCW30MPyiTP1vkrTysnMqJUxho1r4Pl2nkjQu
pUS9oGIH/71Bc4NhLJJ1R29FroQTjkVZnGg830VSSpIuXpF1xR6FrjUNlH7S
bytylpeKO0i5Rz6WRo6SpmN7WF6ZXrCSvCOmmD5JcvZeenZDwhuuQ2b9M2cI
EstEsqY5xwfRpQV8y9ELpD2KRozvsmdQFXhQ27x0iZSzyrNEudAPE8NNe5ka
v0n8ozEFR1UCUPdge79+HA+ejrA7mp2CYMjKezsmx4egw95WbpnBaGnjeZTC
XYPB3fEmyRKhqF3iNkvOCyAjT/y58vHzsKnw6C3VFISVYKlhY04YWqOeV3TZ
CaN6rl2HkI3ABpJSu83iVvi7u1fEe13E990pX9smtkHhLaEhV3gJES9+nz4/
JY6UcpEd9ipX4BJ051OIqjTMpIuKpk5FweRu4BU3iiyutNwWeuEUjasy3Jdv
++Sgraue7JpJRa32+ls6JorIF/1v3DWJWGjq7w/fG70DlzqyWEefxsLSiiYg
PqGKo7SRVYrHaYjRVlkX46CGjYvBUfU5cXjjo4tscw2E92SEz8gm4kBqAPun
pSYckwara58nwq11ViFZSgsapiqLsV7XYDtr1i37XQvjlTXB0oXXoR7JV6lK
14H9H/2+diopCF0d6i3ktA/q1oolkq9ToNkW4TNhvh/8fP6fuadS0a7hV04z
MoBH3fKTJqc0NIp7uArpGkFUNU3DJTMdteXUNb0lxqDyV/4nUSIDjMHySE0Z
M+y0SE+oCMoA5dK3Pf4vzMCmBU5xXmdgW1NxAd3HHEp86mUaYTeRP6kUb1lj
gZZ/ijR62Jng8o8qhAVZYeNCpo5s0i2MiFtkH5A7DObe85mWY/++UK5OKvUw
3IMcl4L330YbhZLfMXYY9FnBOtbkbex54xwpXmYPh8kEQ1+zHNR80D/bauqI
7MGLJszlu6i3VrlfXngUB6Rci4aO5GGjwIMv17c4Y/qY3hr+8jxIUjbxWcc6
Q2Ymgly+Kx9lblSx06ddEHhGV5lnIB8JGlJOawyOQcSb/rBYvwA8yRs1hNa7
w8X1Hn3+axHCw24aRsqTQXCcNhRIFqOfPCQvdpKMOHd+iqNJMaIv23itLySj
khwoWZ1RjDiO8/A3ZW9sDjCiXdK7D67fozT2bgpexfQdW3v5ATzmlfWSHE9B
85HwWin+UhCNKV8N+/wgsEJ1s5Ch45eEaydr8hyyzzi7q5vT0Rf/t9lBkRAo
9ah7CgfQ0SQqojZmEAJSSS+Z84wOxegLsympIagryWTWaDDQusK7LUk/FLOO
Mw3WbxkwPsRPEUg/jyMPs/4M2VSBBc2M4czj4VC2n5tsvl4FtfxPU4sOV3Dh
7T5kf6+DwdT8Ye2QOIhVCQ/x8jZhZJOPSfzOFFSPmrR3Gr5EyoB7rG71H3jp
m7vddYaFsKZd9ROl2VhUMRIsP/saebrj8in5uVj1o0xou/IbtDwMuyr16B50
VzyrHR+3g0RVk7G/XIT2WO23LTs70ZBnIMX0egDZtC8/SubeoSA9W9/0YSm+
LA0Pvd/DhlkBt1foQzmKpuvfbE7phofgrKhHFQObpR8Q4+M9uG30W02WVoGk
M0yffKUPYLeXsmVNq/CXcdn46ahebPC1OW2rV41QTVbVtuaPIC/EiHpuq8Xw
Y9b1Lvl+rNp3qv662lP88u25o3qRg/U65DEBswaoHVi37FM/gJrbLRGz5o2w
KJA/9URsCDKXjx9/Z/wcF6Q9FD65DcOl4dV8rv0LqM8cEei1HUHi3sCd4uSX
2CI4z8/6OQJGUfnP3BvNWN10h66Y8wnhVmq8lCoW/LTr+9Zbj4Ir9fhsw+Ir
lB06md/41yhWpm7WNeu/xuIXSSFu5hiEuerdPGdbwXS04MswHYeaXdYXseY2
RFnS8qq440gM2TGzdsOfkB+0+lhImwBZVI792b4d+/ViMyQ0uWC3vg62aOyA
m0znD95uLviMZSiFkp1YfJM8GkT9jIsbvgVM271D9mr3UAelSZjmOQt58ndh
yxCzYBVrEtxPgt/aXnZhzWWq077gKVD/UTx74gYbPCXu9FGRaeQ5hKp7unbj
3G16ssjLaeixtVyE9XrwTC2BXzjgC7Yqcjukxd5D2OqceszmGTwt9BmSnHuP
gV3vuiMaZtDEO2UsOPQBjLtQGnefRTynyrmD1YsTdF1h6eVZMLLjNmWXf0Q0
7R51oHQOx7LE5Aae9KF8UvrOZhsCayW1Eg886kdCLs1zyI6A55nQLAajH94O
s8qP7QmEz3cZyjH7IdXQwLB1JrCVcciWt7Ef11KPtySdJHBpmS/hTWc/Tuve
+rohgkBlj+eliIV+KMRutVx3n0CPAUXQyoCDZf3ITb0PCTRtuqzVDg76ZsfY
98sIWCoSX+1MOUhzrjxh/oSAX9gh0lFbDlZ2OoTQGgio/SoR8/fgYKAr9S5f
N4E2+8djVVc5qI1bIne/J+DNulpnnMBBhoGHauFHAguyh9U6kjiwLtaoOjhE
ID9Ee+zzTQ6ehne0UacIZFr8uq/0kIPMXTqp9jMELF47ezAfcRA8muMkTxBI
qAuvMn/CwQ4b/+Hn3wnQbgTo+D3lYDVvT1HaEoGYxjCv5UYORpgGvp4/Cfw9
XbLzRjMHDb6Fu/f8QyCDKZ+u0spBlqzA4soKgd/MpYz6dg7+B7ZtB7E=
"]]}, "Charting`Private`Tag$3205#1"]}, {}, {}}, AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {"k x tube radius", "response"}, AxesOrigin -> {0, 0}, DisplayFunction -> Identity, Frame -> {{False, False}, {False, False}}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImagePadding -> All, Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({(Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[Part[#, 1]], (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({(Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[Part[#, 1]], (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[Part[#, 2]]}& )}}, PlotRange -> {{0, 1.}, {0., 0.43437089488619995`}}, PlotRangeClipping -> True, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}],Placed[\!\(TraditionalForm\`FormBox[TemplateBox[{"\<\"Relative error of (Capsule - Tube)\"\>"}, "\<LineLegend\>", DisplayFunction -> \((FormBox[StyleBox[StyleBox[PaneBox[TagBox[GridBox[{{TagBox[GridBox[{{GraphicsBox[{{Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[0]]], PointSize[0.5], Opacity[1. ], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6]], {LineBox[{{0, 10}, {20, 10}}]}}, {Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[0]]], PointSize[0.5], Opacity[1. ], RGBColor[0.368417, 0.506779, 0.709798], AbsoluteThickness[1.6]], {}}}, AspectRatio -> Full, ImageSize -> {20, 10}, PlotRangePadding -> None, ImagePadding -> Automatic, BaselinePosition -> \((Scaled[0.1] -> Baseline)\)], #}}, GridBoxAlignment -> {"\<Columns\>" -> {Center, Left}, "\<Rows\>" -> {{Baseline}}}, AutoDelete -> False, GridBoxDividers -> {"\<Columns\>" -> {{False}}, "\<Rows\>" -> {{False}}}, GridBoxItemSize -> {"\<Columns\>" -> {{All}}, "\<Rows\>" -> {{All}}}, GridBoxSpacings -> {"\<Columns\>" -> {{0.5}}, "\<Rows\>" -> {{0.8}}}], "\<Grid\>"]}}, GridBoxAlignment -> {"\<Columns\>" -> {{Left}}, "\<Rows\>" -> {{Top}}}, AutoDelete -> False, GridBoxItemSize -> {"\<Columns\>" -> {{Automatic}}, "\<Rows\>" -> {{Automatic}}}, GridBoxSpacings -> {"\<Columns\>" -> {{1}}, "\<Rows\>" -> {{0}}}], "\<Grid\>"], Alignment -> Left, AppearanceElements -> None, ImageMargins -> {{5, 5}, {5, 5}}, ImageSizeAction -> "\<ResizeToFit\>"], LineIndent -> 0, StripOnInput -> False], {FontFamily -> "\<Arial\>"}, Background -> Automatic, StripOnInput -> False], TraditionalForm] &)\), Editable -> True, InterpretationFunction :> \((RowBox[{"\<LineLegend\>", "\<[\>", RowBox[{RowBox[{"\<{\>", RowBox[{"\<Directive\>", "\<[\>", RowBox[{RowBox[{"\<Opacity\>", "\<[\>", "\<1.`\>", "\<]\>"}], "\<,\>", InterpretationBox[ButtonBox[TooltipBox[GraphicsBox[{{GrayLevel[0], RectangleBox[{0, 0}]}, {GrayLevel[0], RectangleBox[{1, \(-1\)}]}, {RGBColor[0.368417, 0.506779, 0.709798], RectangleBox[{0, \(-1\)}, {2, 1}]}}, AspectRatio -> 1, Frame -> True, FrameStyle -> RGBColor[0.24561133333333335`, 0.3378526666666667, 0.4731986666666667], FrameTicks -> None, PlotRangePadding -> None, ImageSize -> Dynamic[{Automatic, 1.35  CurrentValue["\<FontCapHeight\>"]/AbsoluteCurrentValue[Magnification]}]], "\<RGBColor[0.368417, 0.506779, 0.709798]\>"], Appearance -> None, BaseStyle -> {}, BaselinePosition -> Baseline, DefaultBaseStyle -> {}, ButtonFunction :> With[{Typeset`box$ = EvaluationBox[]}, If[Not[AbsoluteCurrentValue["\<Deployed\>"]], SelectionMove[Typeset`box$, All, Expression]; FrontEnd`Private`$ColorSelectorInitialAlpha = 1; FrontEnd`Private`$ColorSelectorInitialColor = RGBColor[0.368417, 0.506779, 0.709798]; FrontEnd`Private`$ColorSelectorUseMakeBoxes = True; MathLink`CallFrontEnd[FrontEnd`AttachCell[Typeset`box$, FrontEndResource["\<RGBColorValueSelector\>"], {0, {Left, Bottom}}, {Left, Top}, "\<ClosingActions\>" -> {"\<SelectionDeparture\>", "\<ParentChanged\>", "\<EvaluatorQuit\>"}]]]], BaseStyle -> Inherited, Evaluator -> Automatic, Method -> "\<Preemptive\>"], RGBColor[0.368417, 0.506779, 0.709798], Editable -> False, Selectable -> False], "\<,\>", RowBox[{"\<AbsoluteThickness\>", "\<[\>", "\<1.6`\>", "\<]\>"}]}], "\<]\>"}], "\<}\>"}], "\<,\>", RowBox[{"\<{\>", #, "\<}\>"}], "\<,\>", RowBox[{"\<LegendMarkers\>", "\<->\>", "\<None\>"}], "\<,\>", RowBox[{"\<LabelStyle\>", "\<->\>", RowBox[{"\<{\>", "\<}\>"}]}], "\<,\>", RowBox[{"\<LegendLayout\>", "\<->\>", "\<\"Column\"\>"}]}], "\<]\>"}] &)\)], TraditionalForm]\),After]]


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
