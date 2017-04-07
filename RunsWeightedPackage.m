(* ::Package:: *)

BeginPackage["RunsWeightedPackage`",{"HypothesisTesting`"}]
(* import package as follows: *)
(*
SetDirectory["/home/pcl312/beaujean/workspace/pValue/"]
 Needs["RunsWeightedPackage`"] 
*)

pValueRuns::usage =
        "pValueRuns[\!\(\*SubscriptBox[\"T\", \"obs\"]\),n] computes the p value \!\(\*FormBox[
RowBox[{\"P\", \"(\", 
RowBox[{
RowBox[{\"T\", \">\", SubscriptBox[\"T\", \"obs\"]}], \"|\", \"n\"}], \")\"}],
TraditionalForm]\) analytically. 
pValueRuns automatically threads over a list of \!\(\*SubscriptBox[\"T\", \"obs\"]\).
For negative \!\(\*SubscriptBox[\"T\", \"obs\"]\), pValueRuns returns Null.
Example: pValueRuns[{6.8, 10.4,15.5},5] gives
{0.05074187128132179`,0.010054846546659557`,0.0009877673221849737`}

"

pValueRunsMC::usage = 
        "pValueRunsMC[\!\(\*SubscriptBox[\"T\", \"obs\"]\),n] computes the p value \!\(\*FormBox[
RowBox[{\"P\", \"(\", 
RowBox[{
RowBox[{\"T\", \">\", SubscriptBox[\"T\", \"obs\"]}], \"|\", \"n\"}], \")\"}],
TraditionalForm]\) by simulating datasets in the Monte Carlo (MC) approach. 
pValueRunsMC automatically threads over a list of \!\(\*SubscriptBox[\"T\", \"obs\"]\).
For negative \!\(\*SubscriptBox[\"T\", \"obs\"]\), pValueRuns returns Null.
Example: pValueRunsMC[{6.8, 10.4,15.5},5] may give something similar to
{0.05165669742394585`,0.010357773181166374`,0.0010297920493147483`}
"

runsSuccess::usage =
"Calculate T, the runs  statistic from  a list of values interpreted as 
independent samples from a standard normal distribution.
Example: runsSuccess[{-1,1,3,-2}]\[Equal]10"

MCexperiments::usage = 
"Number of pseudo experiments used to estimate p-value in the Monte Carlo approach. Default: 5000.
Change value using RunsWeightedPackage`MCexperiments = n"

MCseed::usage =
"Random number seed used to generate pseudo experiments. Default: 134143"

ECDF::usage =
"
Computes the empirical cumulative distribution function for a list of samples 
Example: ECDF[{1,1.3,2.4,4.1} ] gives
{{1, 1/4}, {1.3, 1/2}, {2.4, 3/4}, {4.1, 1}}
"

Begin["`Private`"]


(*our model function that gives weight to an individual run of length l, P (L<\[Chi]^2|l) . *)
 Pmodel[x_,runLength_Integer]:=CDF[ChiSquareDistribution[runLength],x] 
DistributeDefinitions[Pmodel] 


(*Probability for one partition \[Pi]=(Subscript[r, 1],... Subscript[r, l , ]... Subscript[r, max]), P (L < \[Chi]^2 | \[Pi])=\!\(
\*SubscriptBox[\(\[Product]\), \(\(l\)\(\ \ \)\)]\(\([\ P \((L < \[Chi]^2 | l)\)\ ]\)^Subscript[r, \ l]\)\)*)
Ppartition[x_, hist_] := Product[ (Pmodel[x, l])^hist[[l]], {l,1,Length[hist]} ];
DistributeDefinitions[Ppartition]


(*create a histogram from the partition *)
HistoFromPartition[part_]:= BinCounts[part,{1,Max[part]+1,1}]
DistributeDefinitions[HistoFromPartition]


(*Multiplicity for one partition, barring exceptions it is the number of permuting the success runs times the number the failures may be shuffled around such that we still have valid runs*)
Wnumber[histo_,m_Integer, r_Integer, n_Integer] :=If[m>n-r+1||r>n,0,(Product[Factorial[histo[[i]]],{i,1,Length[histo]}])^-1*Gamma[n-r+2]/Gamma[n-r+2-m]]
DistributeDefinitions[Wnumber]


(*Probability for all partitions {\[Sigma]}, P (L < \[Chi]^2 | \[Pi])=\!\(
\*SubscriptBox[\(\[Product]\), \(\(l\)\(\ \ \)\)]\(\([\ P \((L < \[Chi]^2 | l)\)\ ]\)^Subscript[r, \ l]\)\)*)
Psigma[x_, m_Integer, r_Integer, n_Integer] :=Module[{ partitions, permMultiplicities, totalMultiplicity},
(* create histogram for every Integer partition *)
partitions = (HistoFromPartition[#])&/@ IntegerPartitions[r,{m}] ;
(* find the multiplicity of each ordered partition *)
permMultiplicities = (Wnumber[#,m,r,n])&/@ partitions;
totalMultiplicity = Total[permMultiplicities];
(* no runs, no contribution *)
If[totalMultiplicity==0, Return[0] ];
(* law of total P. *)
Return[Sum[ permMultiplicities[[i]] *Ppartition[x,partitions[[i]] ] , {i,1,Length[partitions]}] ]
]
DistributeDefinitions[Psigma]


(*Put it all together P (L<\[Chi]^2|N) to get the p - value*)
finalP[x_, n_] :=Module[{R,p},
(*R= Rnumber[n];*)
p=1-(2^n-1)^-1  ParallelSum[Psigma[x,  m,r,n],{r,1,n},{m,1,Min[r,n-r+1] },Method->Automatic];
Return[p]
]
DistributeDefinitions[finalP]


(*Allow for multiple \[Chi]^2 values, compute the function P (L<\[Chi]^2|N) only once, keeping \[Chi]^2 arbitrary. Can't do parallel replacement (why not?), so define a nested function that can be parallelized.*)
pValueRuns[listX_, n_Integer] := Module[{function, returnValues},
(*returnValues=Null;*)
(* check for proper numerical input: only Subscript[T, obs]>= 0 is acceptable *)
If[Length[listX]==0,
If[ (*NumericQ[listX] &&*) FreeQ[listX,x_/;x<  0], Return[finalP[listX,n]],Return[Null] ]
];
If[Length[listX]!=0&& ( (*!VectorQ[listX,NumericQ] ||*) !FreeQ[listX,x_/;x<  0] ),Return[]];

function =finalP[x,n]//Evaluate; 
PofValue[z_]:= Evaluate[function/.{x-> z}];(*reuse stored definition *)
DistributeDefinitions[PofValue];
returnValues = Parallelize[Map[PofValue,listX],Method-> "CoarsestGrained"];

Return[ returnValues]
]/;n>= 2

DistributeDefinitions[pValueRuns]


runsSuccess[listSamples_]:= Module[{temp},

(* check for proper input *)
If[Length[listSamples]==0 || !VectorQ[listSamples,NumericQ] ,Return[Null] ];

(* split into runs *)
temp = Split[listSamples,  #1 * #2 >0 &];
temp =( temp /.{x_ :>  0 /;x<=0})^2; (* remove failures *)
temp = Map[Total,temp]; (* sum in each run *)
Return[temp//Max];
];


(*ECDF *)
ECDF[dataList_]:=Module[{sample, points,nPoints},
sample=Tally[dataList];
nPoints=Length[dataList];
points={Transpose[sample][[1]],Accumulate[Transpose[sample][[2]]]/nPoints};
Return[Transpose[points]]
]


(**** MC part ****)
MCexperiments = 5000;
MCseed = 134143;


pValueRunsMC[runsStatistics_, n_Integer] := Module[
{listX,k,badIndices,fakeExperiments, returnValues, chi2Runs,pointsECDF, ECDFinterpol,pValue},

(* turn single argument into a list, so from now on we only deal with a list of input values *)
If[Length[runsStatistics]==0,listX={runsStatistics}, listX=runsStatistics];
k=Length[listX];


(* check for proper numerical input: only Subscript[T, obs]>= 0 is acceptable *)

(*for input {3}, Position[...] returns {{0},{}} *)
badIndices= Drop[Position[listX,x_/;!NumericQ[x] || x<0]//Flatten, 1 ];
If[ Length[listX]>0 && (Length[badIndices]==Length[listX]),
	Print["Error: all arguments are either negative or not numerical quantities!"];
	Return[];
];

If[  (Length[badIndices]>0) && Length[badIndices]<Length[listX] ,
	Print["Warning: some arguments are either negative or not numerical quantities!"];
];


(* initialization  *)
SetSharedVariable[n]; (* make it available to all kernels *)


(* generate data and extract runs success and failure statistic for each experiment  *)
SeedRandom[MCseed,Method-> "MersenneTwister"];
fakeExperiments = Table[RandomReal[NormalDistribution[0,1],n],{i,1,MCexperiments}]; 

chi2Runs=DeleteCases[ParallelMap[runsSuccess, fakeExperiments, Method-> "CoarsestGrained"],x_/;x==0]//Sort;
Clear[fakeExperiments];

(* compute ECDF *)
pointsECDF = ECDF[chi2Runs ];

(* build interpolation of ECDF*)
ECDFinterpol = Interpolation[pointsECDF,Method-> "Spline", InterpolationOrder->1];

(* a list of values *)
returnValues= Table[ 
If[ NumericQ[listX[[i]]]  ,
	If[listX[[i]]>0, (* now we have valid arguments *)
		pValue =1 - ECDFinterpol[listX[[i]]];
		(* check if value inside polynomial range *)
		If[listX[[i]] < Min[chi2Runs],Print["warning: observed value "<>ToString[listX[[i]],InputForm]<>" at position "<>ToString[i]<>" is smaller than any simulated value. For more accuracy, increase MCexperiments"];pValue=0];
		If[listX[[i]] > Max[chi2Runs],Print["warning: observed value "<>ToString[listX[[i]],InputForm]<>" at position "<>ToString[i]<>" is larger than any simulated value. For more accuracy, increase MCexperiments"];pValue=1];
		pValue,  (*return value of if-true branch*)

		Null(* else invalid *)
	],

	Null
],{i,k}];

Return[returnValues];

]



(*MC test cases *)

(*
(* negative input *)
pValueRunsMC[-1.3,2]
(* non numeric input *)
pValueRunsMC[,2]
pValueRunsMC["a",2]
pValueRunsMC[{"a","b"},2]
pValueRunsMC[{aasd,bfhasf},5]
(* poor ranges *)
pValueRunsMC[{10^-12,6.8, 10.4,"a",1565.5},5]

(* critical values from paper *)
pValueRunsMC[6.8,5]
pValueRunsMC[{6.8, 10.4,15.5},5]
*)


(* Test cases *)
precision=10^-6 

 Pmodel[0,4]==0
Pmodel[10^5.,4]>1-precision

Ppartition[x, {0,0,1,2}]//FullSimplify

(* critical values from paper *)
Abs[finalP[11.5,25]-0.05]<0.001
Abs[finalP[21.6,25]-0.001]<0.001
Abs[finalP[12.8,10]-0.01]<0.001

finalP[x,3]


(* negative input *)
pValueRuns[-1.3,2]
(* non numeric input *)
pValueRuns[{"a","B"},2]
pValueRuns[{"a",1},2]
pValueRuns[{a,b},5]
(* critical values from paper *)
pValueRuns[6.8,5]
pValueRuns[{6.8, 10.4,15.5},5]





(* data splitting test cases *)
runsSuccess[{1,-1,-3,2}]==4
runsSuccess[-{1,-1,-3,2}]==10
runsSuccess[{1,-1,-3,a}]==Null
runsSuccess[{1,2,-1,1,1,1,-2}]==5
runsSuccess[-{1,2,-1,1,1,1,-2,-3}]==13
runsSuccess[{}]==Null



End[ ] (* "`Private`" *)

EndPackage[ ]









