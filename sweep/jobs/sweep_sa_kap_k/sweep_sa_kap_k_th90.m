(* ::Package:: *)

(* ::Input:: *)

LaunchKernels[64];

lmat=Get["cache/lmat_n16_4_16.m"]


f[x_]:=Ceiling[4.047188292450972` +3.127154434408593` x-0.32286996149984337` x^2+0.015522386535739523` x^3]
index[\[Kappa]_]:=With[{t=Table[i,{i,4,16,3}]},Position[t,First@Nearest[t,f[\[Kappa]]]]]


lmat2=Map[Normal[#]/.{b->\[Kappa] dr}&,lmat];


lmatfun[n_][\[CapitalTheta]_,dr_,dt_,\[Sigma]a_,\[Sigma]m_,k_,\[Kappa]_]=lmat2[[n]];


eval[\[CapitalTheta]_,dr_,dt_,\[Sigma]a_,\[Sigma]m_,k_,\[Kappa]_]:=Eigenvalues[N[lmatfun[First@index[\[Kappa]]][\[CapitalTheta],dr,dt,\[Sigma]a,\[Sigma]m,k,\[Kappa]]]]



myparams={
	dr->0.05,
	dt->0.1,
	sm->0.01
}

res90=Parallelize@Table[
{{\[Pi]/2, dr, dt, sa, sm dr \[Kappa], k, \[Kappa]} /. myparams, eval@@({\[Pi]/2, dr, dt, 2 sa, sm dr \[Kappa], k, \[Kappa]} /. myparams)},
{sa,Table[i, {i, -3, 3, 0.05}]},
{\[Kappa],Table[i, {i, 0, 25, 0.1}]},
{k,Join[{0}, Table[1/100 2^i, {i,0,9}]]}
];
Export[$ScriptCommandLine[[3]]<>"/sweep_sa_kap_k_th90_"<>$ScriptCommandLine[[2]]<>".m", res90]
