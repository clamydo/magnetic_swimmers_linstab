(* ::Package:: *)

(* ::Section:: *)
(*Definitions*)


(* ::Subsection:: *)
(*Mathematical primitives*)


(* ::Subsubsection:: *)
(*basis of spherical coordinates.*)


ClearAll[er,e\[Theta],e\[Phi],\[Theta],\[Phi]]
er[\[Theta]_,\[Phi]_]:={Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]}
e\[Theta][\[Theta]_,\[Phi]_]:={Cos[\[Theta]]Cos[\[Phi]],Cos[\[Theta]]Sin[\[Phi]],-Sin[\[Theta]]}
e\[Phi][\[Theta]_,\[Phi]_]:={-Sin[\[Phi]],Cos[\[Phi]],0}


(* ::Subsubsection:: *)
(*gradient operator on unit sphere in spherical coordinates.*)


ClearAll[gradS,divS,diffS,\[Theta],\[Phi],m\[Theta],m\[Phi]]
(*gradS[m\[Theta]_,m\[Phi]_]=(e\[Theta][m\[Theta],m\[Phi]]D[#,\[Theta]]+e\[Phi][m\[Theta],m\[Phi]] 1/Sin[m\[Theta]] D[#,\[Phi]])&;
divS[m\[Theta]_,m\[Phi]_]=(e\[Theta][m\[Theta],m\[Phi]].D[#,\[Theta]]+1/Sin[m\[Theta]] e\[Phi][m\[Theta],m\[Phi]].D[#,\[Phi]])&;
diffS[\[Theta]_,\[Phi]_][f_? MatchQ[#,SphericalHarmonicY[l_,m_,\[Theta]_]&]:=
diffS[\[Theta]_,\[Phi]_]=divS[\[Theta],\[Phi]][gradS[\[Theta],\[Phi]][#]]&;*)
gradS=(e\[Theta][\[Theta],\[Phi]]D[#,\[Theta]]+e\[Phi][\[Theta],\[Phi]] 1/Sin[\[Theta]] D[#,\[Phi]])&;
divS=(e\[Theta][\[Theta],\[Phi]].D[#,\[Theta]]+1/Sin[\[Theta]] e\[Phi][\[Theta],\[Phi]].D[#,\[Phi]])&;
diffS[f_]:=With[{l=First[List@@f]}, -l(l+1) f]/;MatchQ[Head@f,SphericalHarmonicY]
diffS[f_]:=divS[gradS[f]];


(*(*test*)*)
(*diffS[\[Psi][\[Theta],\[Phi]]]==1/Sin[\[Theta]] D[Sin[\[Theta]]D[\[Psi][\[Theta],\[Phi]],\[Theta]],\[Theta]]+1/Sin[\[Theta]]^2 D[\[Psi][\[Theta],\[Phi]],{\[Phi],2}]//FullSimplify*)
(*diffS[SphericalHarmonicY[a+b,m,\[Theta],\[Phi]]]*)


(* ::Subsubsection:: *)
(*Diffusion operator*)


diff[\[Theta]_,\[Phi]_]=(dr With[{r={x,y,z}},Div[Grad[#,r],r]]+ds diffS[\[Theta],\[Phi]][#])&;


(* ::Subsubsection:: *)
(*Gradient operator*)


ClearAll[nablas,nablav,nablat,s]
nablas[s_]:=With[{r={x,y,z}},Grad[s,r]]
nablav[s_]:=With[{r={x,y,z}},Table[D[s[[i]],r[[j]]],{i,1,3},{j,1,3}]]
nablat[s_]:=With[{r={x,y,z}},Table[Div[s[[All,j]],r],{j,1,3}]]


(* ::Subsubsection:: *)
(*Double dot product of matrices*)


ClearAll[ddot]
ddot[a_?MatrixQ,b_?MatrixQ]:=Tr[b.Transpose[a]];


(* ::Subsubsection:: *)
(*Dyadic product*)


ClearAll[dyad]
dyad[a_?VectorQ,b_?VectorQ]:=Outer[Times,a,b]


(* ::Subsubsection:: *)
(*Tangential projector*)


pr=IdentityMatrix[3]-dyad[er[\[Theta],\[Phi]],er[\[Theta],\[Phi]]];


(* ::Section:: *)
(*Spherical harmonics*)


(* ::Subsection:: *)
(*Calculate projection coefficients*)


ClearAll[indeces,coeff,coeffList,combine]
indeces[n_]:=Flatten[Table[{#,m},{m,-#,#}]&/@Table[l,{l,0,n}],1]
coeff[f_+g_,{l_,m_}]:=coeff[f,{l,m}]+coeff[g,{l,m}]
coeff[f_,{l_,m_}]:=Integrate[Sin[\[Theta]] f Conjugate@SphericalHarmonicY[l,m,\[Theta],\[Phi]],{\[Theta],0,\[Pi]},{\[Phi],0,2\[Pi]}]
(*coeff[c_ f_,n_]:=c coeff[f,n]/;(FreeQ[c,\[Theta]]\[And]FreeQ[c,\[Phi]])*)
coeffList[f_,n_]:=MapThread[List,{
	indeces[n],
	Parallelize[Map[
			coeff[f,#]&,
			indeces[n]
		]
	]
}]

combine[c_]:=Total@Map[Last@# SphericalHarmonicY[#[[1,1]],#[[1,2]],\[Theta],\[Phi]]&,Select[c,Last@#=!=0&]]
(*combine[c_]:=Total@Map[Last@# SphericalHarmonicY[#[[1,1]],#[[1,2]],\[Theta],\[Phi]]&,c]*)
combineHold[c_]:=Total@Map[Last@# Hold[SphericalHarmonicY][#[[1,1]],#[[1,2]],\[Theta],\[Phi]]&,Select[c,Last@#=!=0&]]
combineSymbolic[c_]:=Total@Map[Last@# sphericalHarmonicY[#[[1,1]],#[[1,2]],\[Theta],\[Phi]]&,Select[c,Last@#=!=0&]]


(* ::Subsection:: *)
(*Orthogonality*)


ClearAll[scalarProduct]
(*SetAttributes[scalarProduct,HoldAll]*)

scalarProduct[f_+g_,{l_,m_}]:=scalarProduct[f,{l,m}]+scalarProduct[g,{l,m}]
scalarProduct[c_ f_,{l_,m_}]:=c scalarProduct[f,{l,m}]/;FreeQ[c,sphericalHarmonicY]

scalarProduct[f_,{l_,m_}]:=With[{a=List@@f},
	KroneckerDelta[a[[1]],l] KroneckerDelta[a[[2]],m]
]/;MatchQ[Head@f,sphericalHarmonicY]

scalarProduct[Power[f_,2],{l_,m_}]:=With[{a=List@@f},
	(-1)^m Sqrt[((2 a[[1]]+1)(2 a[[1]]+1)(2 l+1))/(4 \[Pi])] Hold[ThreeJSymbol[{a[[1]],0},{a[[1]],0},{l,0}] ThreeJSymbol[a[[1;;2]],a[[1;;2]],{l,-m}]]
]/;MatchQ[Head@f,sphericalHarmonicY]\[And]MatchQ[Head@f,sphericalHarmonicY]

scalarProduct[f_ g_,{l_,m_}]:=With[{a=List@@f,b=List@@g},
	(-1)^m Sqrt[((2 a[[1]]+1)(2 b[[1]]+1)(2 l+1))/(4 \[Pi])] Hold[ThreeJSymbol[{a[[1]],0},{b[[1]],0},{l,0}] ThreeJSymbol[a[[1;;2]],b[[1;;2]],{l,-m}]]
]/;MatchQ[Head@f,sphericalHarmonicY]\[And]MatchQ[Head@g,sphericalHarmonicY]

(*scalarProduct[f_,{l_,m_}]:=With[{a=List@@f},
	2 Sqrt[\[Pi]] KroneckerDelta[0,l] KroneckerDelta[0,m]
]/;FreeQ[Head@f,sphericalHarmonicY]*)


(* ::Section:: *)
(*L coefficients helper*)


(* ::Text:: *)
(*Define a function to map the four dimensional tensor Subscript[L^(m,n), l,j] to a two dimensional tensor Subscript[L, r(j,n),c(l,m)].*)


ClearAll[mode2position,position2mode]
mode2position[{l2_,m2_,l3_,m3_}]:=If[l2>=0\[And]-l2<=m2<=l2\[And]l3>=0\[And]-l3<=m3<=l3,
    {l3^2+(m3+l3)+1,l2^2+(m2+l2)+1}]
position2mode[{x_,y_}]:=Block[{l2,m2,l3,m3},
    l2=Floor@Sqrt[y-1];
    m2=y-l2^2-l2-1;
    l3=Floor@Sqrt[x-1];
    m3=x-l3^2-l3-1;
    {l2,m2,l3,m3}
]


(* ::Input:: *)
(*(*Check inverse function*)*)
(*With[{n=1},*)
(*Flatten[*)
(*	Table[*)
(*		Flatten[*)
(*		Table[{l2,m2,l3,m3},{l2,0,n},{m2,-l2,l2}],1],*)
(*{l3,0,n},{m3,-l3,l3}],*)
(*1]*)
(*==Table[position2mode[{l,m}],{l,1,(n+1)^2},{m,1,(n+1)^2}]*)
(*]*)


(* ::Input:: *)
(*With[{n=16},*)
(*	Flatten@Table[mode2position[position2mode[{l,m}]]=={l,m},{l,1,n},{m,1,n}]*)
(*]//DeleteDuplicates*)
