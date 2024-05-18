(* ::Package:: *)

<<"definitions.m"


(* ::Code:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Section:: *)
(*Variables*)


(* ::Text:: *)
(*Assume the perturbation to be in the x-z-plane.*)


vk=k er[\[CapitalTheta],0]


(* ::Subsubsection:: *)
(*Magnetic field*)


(* ::Text:: *)
(*Assume the magnetic field to be aligned with the z-axis.*)


ClearAll[b,vb]
vb={0,0,b}


(* ::Subsubsection:: *)
(*Fixpoint*)


ClearAll[\[Psi]0]
\[Psi]0[\[Theta]_,\[Phi]_]:=\[Kappa]/(4\[Pi] Sinh[\[Kappa]]) Exp[\[Kappa] er[\[Theta],\[Phi]].{0,0,1}]


\[Psi]0[\[Theta],\[Phi]]


(* ::Subsubsection:: *)
(*Transformed general stress tensor*)


ClearAll[\[CapitalSigma]sym]
\[CapitalSigma]sym=Map[Indexed[\[Sigma],#]&,Table[{i,j},{i,1,3},{j,1,3}],{2}]


(* ::Subsubsection:: *)
(*Transformed of flow-field*)


ClearAll[u]
u=I/k^2 (IdentityMatrix[3]-1/k^2 dyad[vk,vk]).\[CapitalSigma]sym.vk


(* ::Subsubsection:: *)
(*Transformed strain and vorticity tensor*)


ClearAll[w,e]
w=I/2 (dyad[vk,u]-dyad[u,vk])//FullSimplify;
%//MatrixForm
e=I/2 (dyad[vk,u]+dyad[u,vk])//FullSimplify;
%//MatrixForm


(* ::Section:: *)
(*Calculate terms of linear integro-differential-operator L*)


(* ::Text:: *)
(*L[Subscript[Y, Subscript[l, 2]]^Subscript[m, 2]]=(-i n\[CenterDot]k + 2 n\[CenterDot]B - (\[DoubleStruckOne]-nn)\[CenterDot]B\[CenterDot]Subscript[\[Del], S] + Subscript[d, S] Subscript[\[CapitalDelta], S] - Subscript[d, t]k^2) Subscript[Y, Subscript[l, 2]]^Subscript[m, 2] + 3 \[Gamma] nn:E'[Subscript[Y, Subscript[l, 2]]^Subscript[m, 2]] Subscript[\[Psi], 0] - (\[DoubleStruckOne]-nn)\[CenterDot][(\[Gamma] E' - W')[Subscript[Y, Subscript[l, 2]]^Subscript[m, 2]]\[CenterDot]n]\[CenterDot]Subscript[\[Del], S]Subscript[\[Psi], 0]*)


terms=<<"cache/result_terms.m";


(* ::Subsection:: *)
(*Calculation*)


(* ::Input:: *)
(*terms=Table[Null,{7}]*)


ClearAll[\[Psi]]
\[Psi]=SphericalHarmonicY;


(* ::Subsubsection:: *)
(*-i n\[CenterDot]k*)


(* ::Input:: *)
(*terms[[1]]=-I er[\[Theta],\[Phi]].vk \[Psi][l,m,\[Theta],\[Phi]];*)
(*%//TraditionalForm*)


(* ::Subsubsection:: *)
(*2 n \[CenterDot] B \[Psi]'*)


(* ::Input:: *)
(*terms[[2]]=2 er[\[Theta],\[Phi]].vb \[Psi][l,m,\[Theta],\[Phi]]//FullSimplify;*)
(*%//TraditionalForm*)


(* ::Subsubsection:: *)
(*- (1-nn) \[CenterDot] B \[CenterDot] Subscript[\[Del], S]\[Psi]'*)


(* ::Input:: *)
(*FullSimplify[-(pr.vb).gradS[\[Psi][l,m,\[Theta],\[Phi]]],m\[Element]Integers\[And]l\[Element]Integers\[And]-l<=m<=l\[And]l>=0]*)


(* ::Input:: *)
(*terms[[3]]=b m Cos[\[Theta]] SphericalHarmonicY[l,m,\[Theta],\[Phi]];*)
(*terms[[4]]=b E^(-I \[Phi]) \[Sqrt]((1+l+m) (l-m)) Sin[\[Theta]] SphericalHarmonicY[l,1+m,\[Theta],\[Phi]];*)
(*terms[[3]]+terms[[4]]//TraditionalForm*)


(* ::Subsubsection:: *)
(*Subscript[d, S] Subscript[\[CapitalDelta], S]\[Psi]' - dt k^2\[Psi]'*)


(* ::Input:: *)
(*terms[[5]]=dr diffS[\[Psi][l,m,\[Theta],\[Phi]]]-dt k^2 \[Psi][l,m,\[Theta],\[Phi]]//FullSimplify;*)
(*%//TraditionalForm*)


(* ::Subsubsection:: *)
(* 1/Subscript[\[Psi], 0] * -(1-nn) \[CenterDot] [(\[Gamma]E' - W') \[CenterDot] n] \[CenterDot] Subscript[\[Del], S]Subscript[\[Psi], 0]*)


(* ::Input:: *)
(*terms[[6]]=-(pr.((\[Gamma] e-w).er[\[Theta],\[Phi]])).FullSimplify[gradS[\[Psi]0[\[Theta],\[Phi]]]/\[Psi]0[\[Theta],\[Phi]]];*)
(*%//TraditionalForm*)


(* ::Subsubsection:: *)
(*1/Subscript[\[Psi], 0] * 3\[Gamma] nn:E' Subscript[\[Psi], 0]*)


(* ::Input:: *)
(*terms[[7]]=3 \[Gamma] ddot[dyad[er[\[Theta],\[Phi]],er[\[Theta],\[Phi]]],e];*)
(*%//TraditionalForm*)


(* ::Subsubsection::Closed:: *)
(*Save terms*)


(* ::Input:: *)
(*terms>>"cache/result_terms.m"*)


(* ::Section:: *)
(*All of the terms, except Subscript[\[Psi], 0], can be expressed with spherical harmonics up to l=3*)


(* ::Subsection:: *)
(*Load coefficients for spherical harmonics expansion*)


shCoeffs=<<"cache/result_shCoeffs.m";


(* ::Subsection:: *)
(*Calculate coefficients for spherical harmonics expansion*)


(* ::Code:: *)
(*LaunchKernels[8];*)


(* ::Input:: *)
(*ClearAll[shCoeffs]*)
(*shCoeffs=With[{expansionOrder={1,1,1,1,0,3,3}},*)
(*MapThread[FullSimplify[coeffList[Select[#1,FreeQ[#,SphericalHarmonicY]&]//ExpandAll//TrigExpand,#2]]&,{terms,expansionOrder}]*)
(*];*)


(* ::Code:: *)
(*shCoeffs*)


(* ::Subsubsection:: *)
(*Make consistency check by comparing expansion with original:*)


(* ::Input:: *)
(*MapThread[combine[#1/.{\[CapitalTheta]->0}]==Select[#2/.{\[CapitalTheta]->0},FreeQ[#,SphericalHarmonicY]&]&,{shCoeffs,terms}]//FullSimplify*)


(* ::Input:: *)
(*MapThread[combine[#1]==Select[#2,FreeQ[#,SphericalHarmonicY]&]&,{shCoeffs,terms}]//FullSimplify*)


(* ::Subsubsection:: *)
(*Save coefficients for spherical harmonics expansion*)


(* ::Input:: *)
(*shCoeffs>>"cache/result_shCoeffs.m"*)


(* ::Section:: *)
(*Calculate truncated expansion of fixpoint*)


psi0expansion=<<"cache/result_psi0expansion_50.m";


(* ::Subsection:: *)
(*Calculation*)


(* ::Input:: *)
(*(*calculate*)*)
(*psi0expansion=Parallelize@Table[{{l,0},coeff[\[Kappa]/(4\[Pi] Sinh[\[Kappa]]) Exp[\[Kappa] Cos[\[Theta]]],{l,0}]},{l,0,50}];*)


(* ::Input:: *)
(*(*save*)*)
(*psi0expansion>>"cache/result_psi0expansion_50.m";*)


(* ::Section:: *)
(*Calculate stress tensor*)


\[CapitalSigma]HDkern=dyad[er[\[Theta],\[Phi]],er[\[Theta],\[Phi]]]-1/3 IdentityMatrix[3]//FullSimplify;


\[CapitalSigma]HDkernCoeffs=Map[coeffList[#,2]&,\[CapitalSigma]HDkern,{2}];


\[CapitalSigma]HDkernSH=Map[combineSymbolic,\[CapitalSigma]HDkernCoeffs,{2}];


(* ::Input:: *)
(*(*test*)*)
(*Map[combine,\[CapitalSigma]HDkernCoeffs,{2}]==\[CapitalSigma]HDkern//FullSimplify*)


\[CapitalSigma]Bkern=FullSimplify[dyad[er[\[Theta],\[Phi]],Normalize[vb]]-dyad[Normalize[vb],er[\[Theta],\[Phi]]],Assumptions->b>0];


\[CapitalSigma]BkernCoeffs=Map[coeffList[#,1]&,\[CapitalSigma]Bkern,{2}];


\[CapitalSigma]BkernSH=Map[combineSymbolic,\[CapitalSigma]BkernCoeffs,{2}];


(* ::Input:: *)
(*(*test*)*)
(*Map[combine,\[CapitalSigma]BkernCoeffs,{2}]==\[CapitalSigma]Bkern//FullSimplify*)


(* ::Text:: *)
(*The calculation of the angular expectation value makes use of the scalar product. For this, the base vectors of the propability distribution function have to be complex conjugated first: Conjugate[SphericalHarmonicY[l,m,\[Theta],\[Phi]]]==(-1)^mSphericalHarmonicY[l,-m,\[Theta],\[Phi]]. But since \[CapitalSigma] is real, \[CapitalSigma]=\[CapitalSigma]**)


ClearAll[\[CapitalSigma]]
\[CapitalSigma]=With[{stress=\[Sigma]a/2 \[CapitalSigma]HDkernSH+\[Sigma]m/2 \[CapitalSigma]BkernSH//FullSimplify},
Map[scalarProduct[#,{l,m}]&,stress,{2}]
];


(* ::Section:: *)
(*Calculate coefficients of linear operator*)


(* ::Subsection:: *)
(*Calculate lCoeff*)


(* ::Text:: *)
(*Minimal l-order is 3, because terms[[1;;5]] need up to order 3*)


coefficientMatrix[order_,\[Psi]order_,subs_] := Module[{\[Psi]0approx,coeff,lCoeffExpr,stresssubs},
	Assert[order>=\[Psi]order]; 
    \[Psi]0approx=combineSymbolic[Map[{First@#,Last@#}&,psi0expansion[[1;;\[Psi]order+1]]]];

	stresssubs=Flatten@MapThread[#1->#2&,{\[CapitalSigma]sym,\[CapitalSigma]},2];
	
    coeff=With[{
	strippedTerms={
			sphericalHarmonicY [l,m,\[Theta],\[Phi]],
			sphericalHarmonicY [l,m,\[Theta],\[Phi]],
			sphericalHarmonicY [l,m,\[Theta],\[Phi]],
			sphericalHarmonicY [l,m+1,\[Theta],\[Phi]],
			sphericalHarmonicY [l,m,\[Theta],\[Phi]],
			\[Psi]0approx,
			\[Psi]0approx
	}},
		MapThread[
			scalarProduct[
				Expand[#1 combineSymbolic[#2],sphericalHarmonicY], 
				{l3,m3}
			]/.stresssubs &,
			{strippedTerms,shCoeffs}
		]
	];
	
	lCoeffExpr = SparseArray@With[{n=order},
					Flatten@Parallelize@Table[
						mode2position[{il2,im2,il3,im3}]->Quiet[ReleaseHold[(Total@coeff/.subs/.{l->il2,m->im2,l3->il3,m3->im3})],{ClebschGordan::phy,ClebschGordan::tri}],
						{il2,0,n},{im2,-il2,il2},{il3,0,n},{im3,-il3,il3}
					]
	];

	lCoeffExpr
];
