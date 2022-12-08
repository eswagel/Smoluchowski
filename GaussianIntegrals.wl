(* ::Package:: *)

BeginPackage["GaussianIntegrals`"];


GaussianIntegral::usage="Integrate from -\[Infinity] to \[Infinity], speeding up Gaussian integrals. Defaults to Integrate when it cannot identify a Gaussian integral.";


Begin["`Private`"]


IsElExp[el_] :=
	(Head[el] == Power && Length[el] > 1 && el[[1]] === E);

IsExp[list_] :=
	IsElExp /@ list

GaussianIntegralResult[a_, b_, c_, n_] :=
	1 / 2 (-a) ^ (-n / 2) E^c (1 / a (-1 + (-1) ^ n) b Gamma[1 + n / 2] 
		Hypergeometric1F1[1 + n / 2, 3 / 2, -(b^2 / (4 a))] + ((1 + (-1) ^ n)
		 Gamma[(1 + n) / 2] Hypergeometric1F1[(1 + n) / 2, 1 / 2, -(b^2 / (4 
		a))]) / Sqrt[-a]);

SingleGaussianIntegral[expr_, intvar_] :=
	Module[{factorlist, exp, coeffs, mask, prefactor, a, b, c, n, const},
		factorlist = Replace[expr, Times -> List, 1, Heads -> True];
		If[Not[Head[factorlist] === List],
			factorlist = {factorlist}
		];
		mask = IsExp[factorlist];
		If[Not[AnyTrue[mask, (# == True)&]],
			(
				Message[SingleGaussianIntegral::mathematica, expr];
				Return[Integrate[expr, {intvar, -Infinity, Infinity}], Module]
			)
		];
		exp = Expand[(Times @@ Pick[factorlist, mask])[[2]]];
		prefactor = Times @@ Pick[factorlist, Not /@ mask];
		coeffs = CoefficientList[exp, intvar];
		If[Length[coeffs] > 3,
			(
				Message[SingleGaussianIntegral::mathematica, expr];
				Return[Integrate[expr, {intvar, -Infinity, Infinity}], Module]
			)
		];
		{c, b, a} = coeffs;
		n = Exponent[prefactor, intvar];
		const = Coefficient[prefactor, intvar, n];
		const * GaussianIntegralResult[a, b, c, n]
	]

SingleGaussianIntegral::mathematica = "Evaluating `1` using Mathematica. Results may be slow.";

GaussianIntegral[expr_, intvar_, simplify_:False] :=
	Module[{termslist, sum},
		termslist = Replace[Expand[expr], Plus -> List, 1, Heads -> True];
		If[Not[Head[termslist] === List],
			termslist = {termslist}
		];
		sum=0;
		Do[sum+=SingleGaussianIntegral[termslist[[i]],intvar];,{i, 1, Length
	[termslist]}];
		
		(*sum = Total[Parallelize[Map[SingleGaussianIntegral[#, intvar]&, termslist
			
			]]];*)
		If[simplify,
			Simplify[sum]
			,
			sum
		]
	]


End[];
EndPackage[];



