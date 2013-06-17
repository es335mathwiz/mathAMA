(* Mathematica Package *)

(* Created by the Wolfram Workbench Jun 17, 2013 *)


(* Mathematica Package *)
BeginPackage["Perturbation`", { "JLink`", "NumericAMA`", "SymbolicAMA`","FunctionCallTree`","Combinatorica`"}]
(* Exported symbols added here with SymbolName::usage *)  

prepAMAModel::usage=
"prepAMAModel[theModel_Symbol,eqns_List,initSSSubs___?OptionQ]"


  xformAMAFO::usage=
"  xformAMAFO[bMat_?MatrixQ,szThetaMat_?MatrixQ,nzCols_List,leadsNeeded_List]"


newAMAModel::usage=
"newAMAModel[model_Symbol,eqns_List]"
AMALags::usage=
"AMALags[model_]"
AMALeads::usage=
"AMALeads[model_]"
AMANeq::usage=
"AMANeq[model_]"
AMAEquations::usage=
"AMAEquations[model_]"


getBDerivs::usage=
"getBDerivs[theModel]"
getLeadsNeeded::usage=
"getLeadsNeeded[theModel]"
getNLags::usage=
"getNLags[theModel]"
getEtaM::usage=
"getEtaM[theModel]"
getEta0::usage=
"getEta0[theModel]"
getEtaP::usage=
"getEtaP[theModel]"
getEta::usage=
"getEta[theModel]"
getNzCols::usage=
"getNzCols[theModel]"
getEpsVars::usage=
"getEpsVars[theModel]"
getNumErrs::usage=
"getNumErrs[theModel]"
getFDerivs02::usage=
"getFDerivs02[theModel]"
getFDerivs03::usage=
"getFDerivs03[theModel]"
isAMAModel::usage=
"isAMAModel[theModel]"
getSplit02::usage=
"getSplit02[theModel,assumedMomSubs]";
getSplit03::usage=
"getSplit03[theModel,assumedMomSubs]";
getSBDerivs::usage=
"getSBDerivs[theModel,assumedMomSubs]"
AMASSGuess::usage=
"AMASSGuess[model]"
ssSolnSubs::usage=
"ssSolnSubs[model,opts]"
computeSteadyState::usage=
"computeSteadyState[model,opts]"
getNewB::usage=
"getNewB[theModel]";
getStochConBs::usage=
"getStochConBs[theModel]"
initSSGuess::usage="steady state gusss options for AMASS"

constructFDerivs::usage=
"constructFDerivs[theModel_Symbol,opts___?OptionQ]"

theIndex::usage="theIndex[beta_List]"


gMults02Java::usage=
"gMults02Java[gPyr01_?MatrixQ,assumedMomSubs_List:{}]"

gMults02Java::usage=
"gMults02Java[gPyr01_?MatrixQ,assumedMomSubs_List:{}]"

gMults03::usage=
"gMults03[gPyr01_?MatrixQ,gArgs_Integer,assumedMomSubs_List:{}]"

notGMults02::usage=
"notGMults02[gPyr01_?MatrixQ,assumedMomSubs_List:{}]"

eqStochSysCon01::usage=
"eqStochSysCon01[theModel_?isAMAModel,assumedMomSubs_List:{}]"

AMAHmat::usage=
"AMAHmat[model_,opts___?OptionQ]"




AMASoln::usage=
"AMASoln[model,opts___] apply numericAMA"


AMACmat::usage=
"AMACmat[model,opts___?OptionQ]"

AMAHmat::usage=
"AMAHmat[model,opts___?OptionQ]"

AMATheta::usage=
"AMATheta[model,opts___?OptionQ]"


infoLag::usage=
"option for AMASoln"

AMABmat::usage=
"AMABmat[model_,opts___?OptionQ]"


AMAQmat::usage=
"AMAQmat[model_,opts___?OptionQ]"

AMASmat::usage=
"AMASmat[model_,opts___?OptionQ]"

AMAS0Inv::usage=
"AMAS0Inv[model_,opts___?OptionQ]"

AMAEigSystem::usage=
"AMAEigSystem[model_,opts___?OptionQ]"



cmpStochLinSys02::usage=
"cmpStochLinSys02[theModel_Symbol,augmentedSubs_List]"

specificKronSum::usage=
"specificKronSum[spRows_List,theLead_Integer,multVals_List,bPows_List]"





eqStochSysCon02::usage=
"eqStochSysCon02[theModel_?isAMAModel,assumedMomSubs_List:{}]"

cmpStochLinSys03::usage=
"cmpStochLinSys03[theModel_Symbol,augmentedSubs_List]"

eqStochSysCon03::usage=
"eqStochSysCon03[theModel_?isAMAModel,assumedMomSubs_List:{}]"


kron::usage=
"kron[a_,b_]:=BlockMatrix[Outer[Times,a,b"



mace::usage=
"mace[exprList_List,limVal_?NumberQ]"
snapShot::usage=
"snapShot[expVal_String]"


(*
If[$OperatingSystem=="Windows",
  ReinstallJava[ClassPath->
"c:\\Documents and Settings\\"<>
 "frbuser\\fromCvs\\stupid\\convergence\\"<>
"src\\mathematica\\pertAMA\\"],
          ReinstallJava[ClassPath->"./"]];*)
If[$OperatingSystem=="Windows",
	AddToClassPath["s:/tryBenchWindows3.5/Perturbation/javaOutput/"],
AddToClassPath["/msu/home/m1gsa00/scratch/tryBenchWindows3.5/Perturbation/javaOutput/"]];

  LoadJavaClass["gov.frb.ma.msu.MultiIndex",StaticsVisible->True,AllowShortContext->True];

defaultOpenFile[]:=$openFiles[[-1]];
resetFile[]:=Map[OpenWrite[Close[#]]& , $openFiles]

closeFile[]:=Map[Close , $openFiles]
fileName[ostr_]:=ostr[[1]]
If[Length[$openFiles]<1,
$openFiles={OpenAppend[
If[StringMatchQ[$Version,"*Windows*"],"c:/temp/defaultTimingFile",
"/tmp/"<>ToString[RunThrough["whoami",""]]<>"defaultTimingFile"]]}];
createNewFile[]:=
With[{newfy=
OpenAppend["timingFile"<>ToString[Abs[Random[]]]<>ToString[SequenceForm[Date[]]]]},
AppendTo[$openFiles,newfy];newfy]
currentFile[]:=defaultOpenFile[]
Global`$globalIndent=""



Begin["`Private`"] (* Begin Private Context *) 








constructNz[incList_List,neq_Integer,nlag_Integer]:=
With[{appearing=Select[incList,#[[2]]<=0&]},
Map[#[[1]]+neq*(nlag+#[[2]]-1)&,appearing]]



constructLeads[incList_List]:=
With[{appearing=Select[incList,#[[2]]>1&]},
Map[{#[[1]],#[[2]]-1}&,appearing]]

constructEta[incList_List,cut_Integer]:=
{Count[incList,{yy_Integer,xx_Integer}/;Or[xx<1,yy>cut]],
Count[incList,{yy_,xx_Integer}/;And[xx==1,yy<=cut]],
Count[incList,{_,xx_Integer}/;xx>1]}

cmpKeepers[leadsNeeded_List,neq_Integer]:=
Union[Flatten[Map[aVarKeep[neq,#]&,leadsNeeded]]]


aVarKeep[neq_Integer,{varNo_Integer,varLeadsNeeded__}]:=
Map[varNo+neq*(#-1)&,{varLeadsNeeded}]






variablesAppearing[incList_List,cut_Integer]:=
With[{justThese=Sort[Union[Flatten[incList,1]],
Or[#1[[2]]<#2[[2]],
And[#1[[2]]==#2[[2]],#1[[2]]=!=1,#1[[1]]<#2[[1]]],
And[#1[[2]]==#2[[2]],#1[[2]]==1,
Or[And[#1[[1]]<=cut,#2[[1]]<=cut],And[#1[[1]]>cut,#2[[1]]>cut]],
#1[[1]]<#2[[1]]],
And[#1[[2]]==#2[[2]],#1[[2]]==1,
And[#1[[1]]>cut,#2[[1]]<=cut]]]&
]},
justThese]



getLag[llPair_List]:=-llPair[[1]]
getLead[llPair_List]:=llPair[[-1]]

stateVariables[modelEquations_]:=
Module[{$AMAjustStateVariablesTime},
Global`AMAPrint[Global`$globalIndent,"timing AMAjustStateVariables"];
(Global`$globalIndent=Global`$globalIndent<>" ");$AMAjustStateVariablesTime=TimeUsed[];
With[{finalResultForTiming=
   Union[Cases[{modelEquations},x_[Global`t]|x_[Global`t+i_]->
   x,Infinity]]},
$AMAjustStateVariablesTime=TimeUsed[]-$AMAjustStateVariablesTime;
Global`AMAPrint[(Global`$globalIndent=StringDrop[Global`$globalIndent,1]),
"AMAjustStateVariables,","done"];
WriteString[defaultOpenFile[],"{AMAjustStateVariables,",
$AMAjustStateVariablesTime,",",MemoryInUse[],",",MaxMemoryUsed[],"}\n"];
finalResultForTiming]]



prepAMAModel[theModel_Symbol,eqns_List,initSSSubs___?OptionQ]:=
With[{epsVars=Union[Cases[eqns,Global`eps[x_][_]->Global`eps[x],Infinity]],
ignore=If[initSSSubs==={},
newAMAModel[theModel,eqns],
newAMAModel[theModel,eqns,initSSSubs]]},
With[{theBMat=AMABmat[theModel],
theS0=-AMAS0Inv[theModel]. AMATheta[theModel]},
With[{fDerivsStuff=constructFDerivs[theModel]},
With[{
fDerivVals=fDerivsStuff[[1]],
nzVals=fDerivsStuff[[2]],
leadsVals=fDerivsStuff[[3]],
etaVals=fDerivsStuff[[4]],
nlagVals=fDerivsStuff[[5]]},
With[{xfrmd=xformAMAFO[theBMat,theS0,nzVals,leadsVals]},
newAMAModel[theModel,fDerivVals,xfrmd[[2]],leadsVals,nlagVals,
etaVals,nzVals,ToExpression[Map["epsXXX" <>ToString[#[[1]]]&,epsVars]]]]]]]]

  xformAMAFO[bMat_?MatrixQ,szThetaMat_?MatrixQ,nzCols_List,leadsNeeded_List]:=
  With[{neq=numRows[bMat],nerrs=numCols[szThetaMat],
    maxLead=Max[Map[Last,leadsNeeded]]},
  With[{nlags=numCols[bMat]/neq}, With[{pth0=If[nlags==1,bMat,
    blockMatrix[{{Drop[IdentityMatrix[neq*nlags],neq]},{bMat}}]],
    keepers=neq*(nlags-1)+Join[Range[neq],neq+cmpKeepers[leadsNeeded,neq]]},
  With[{pth=multExtndMult[bMat,pth0,maxLead]},
  With[{ePth0=If[nlags==1,szThetaMat,
    blockMatrix[{{zeroMatrix[neq*(nlags-1),nerrs]},{szThetaMat}}]]},
  With[{ePth=multExtndMult[bMat,ePth0,maxLead]},
  With[{flatB=blockMatrix[{{pth[[keepers,nzCols]],
        ePth[[keepers,All]]}}]},{
    blockMatrix[{{Drop[pth,neq*(nlags-1)],Drop[ePth,neq*(nlags-1)]}}],
    blockMatrix[{{IdentityMatrix[numCols[flatB]]},{flatB}}]}]]]]]]]



newAMAModel[model_Symbol,eqns_List]:=
  Block[{},Clear[model];AMAEquations[model]^=eqns;isAMAModel[model]^=True;]
newAMAModel[model_Symbol,eqns_List,initSSRule_?OptionQ]:=
  Block[{},Clear[model];
AMAEquations[model]^=eqns;
model/:AMASSGuess[model]=initSSRule;
isAMAModel[model]^=True;
]


newAMAModel[theModel_Symbol,
sFDerivs_List,
bDerivs_?MatrixQ,
leadsNeeded_List,
lags_Integer,
theEta:{etaM_Integer,eta0_Integer,etaP_Integer},
nzCols_List,
epsVars_List]:=
Module[{},
getBDerivs[theModel]^=bDerivs; 
getLeadsNeeded[theModel]^=leadsNeeded;
getNLags[theModel]^=lags;
getEtaM[theModel]^=etaM;
getEta0[theModel]^=eta0;
getEtaP[theModel]^=etaP;
getEta[theModel]^={etaM,eta0,etaP};
getNzCols[theModel]^=nzCols;
getEpsVars[theModel]^=epsVars;
getNumErrs[theModel]^=Length[epsVars];
getFDerivs02[theModel]^=splitFDerivs02[sFDerivs,Apply[Plus,theEta]];
getFDerivs03[theModel]^=splitFDerivs03[sFDerivs,Apply[Plus,theEta]];
isAMAModel[theModel]^=True;
]

isAMAModel[___]:=False;


AMALags[model_]:=(AMALags[model]^=(-1)*lagsLeads[
AMAEquations[model]][[1]])/;AMAEquations[model]=!=AMAundefinedModel;
AMALags[___]:=AMAundefinedModel;

AMANeq[model_]:=Length[AMAEquations[model]]/;
AMAEquations[model]=!=AMAundefinedModel;
AMANeq[___]:=AMAundefinedModel;


AMALeads[model_]:=(AMALeads[model]^=lagsLeads[
AMAEquations[model]][[-1]])/;AMAEquations[model]=!=AMAundefinedModel;
AMALeads[___]:=AMAundefinedModel;


makeSomeMomSubs[epsVars_List]:=
With[{theAlts=Apply[Alternatives,epsVars]},
{(ppp___*xx1_[yy___]*xx2_[zz___]/;
And[MatchQ[xx1,theAlts],MatchQ[xx2,theAlts],xx1=!=xx2]->0.0),
(xx_[___]^zz_Integer/;MatchQ[xx,theAlts]->Global`mom[xx,zz]),
(xx_[yy___]/;MatchQ[xx,theAlts]->Global`mom[xx,1])}]



splitFDerivs02[fDeriv_?MatrixQ,argDim_Integer]:=
With[{nif=numInFacet[argDim,2]},
{fDeriv[[All,Range[argDim]]],fDeriv[[All,argDim+Range[nif]]]}]


splitFDerivs03[fDeriv_?MatrixQ,argDim_Integer]:=
With[{nif=numInFacet[argDim,2]},
{fDeriv[[All,Range[argDim]]],fDeriv[[All,argDim+Range[nif]]]}]

lagsLeads[modelEquations_]:=
      Union[
Join[Cases[modelEquations,x_[Global`t_]->0,Infinity],
Cases[modelEquations,x_[Global`t+v_]->v,Infinity]]];

completeArgList[modelEquations_]:=
Module[{$amacompleteArgListTime},
Global`AMAPrint[Global`$globalIndent,"timing AMAcompleteArgList"];
(Global`$globalIndent=Global`$globalIndent<>" ");$amacompleteArgListTime=TimeUsed[];
With[{finalResultForTiming=
With[{allv=stateVariables[modelEquations],
ll=lagsLeads[modelEquations]},
Flatten[Map[Through[allv[Global`t+#]]&,Range[-getLag[ll],getLead[ll]]]]]},
$amacompleteArgListTime=TimeUsed[]-$amacompleteArgListTime;
Global`AMAPrint[(Global`$globalIndent=StringDrop[Global`$globalIndent,1]),
"AMAcompleteArgList,","done"];WriteString[defaultOpenFile[],
"{AMAcompleteArgList,",$amacompleteArgListTime,",",
MemoryInUse[],",",MaxMemoryUsed[],"}\n"];
finalResultForTiming]]

incompleteArgList[modelEquations_]:=
                Union[Cases[{modelEquations},x_[Global`t]|x_[Global`t+i_],Infinity]]




errorVarCutoff[aModel_Symbol]:=
Module[{},
With[{sv=
stateVariables[AMAEquations[aModel]]},
With[{pos=Flatten[Position[sv,Global`eps[_]]]},
Min[pos]-1]]]



variableExplode[mapper_Function,dim_Integer]:=
With[{allMi=Join[miInFacet[dim,1],miInFacet[dim,2]]},
Map[theIndex,Map[mapper,allMi]]]

makeMatlabFuncs[model_Symbol,opts___?OptionQ]:=
(With[{eqns=AMAEquations[model]},
With[{linPtSubs=Flatten[hmatSSSolnSubs[model,opts]]},
With[{allv=completeArgList[eqns],
sv=stateVariables[eqns]},
With[{varsForSubbing=Table[Unique["subVar"],{Length[allv]}],
appears=Map[incompleteArgList,eqns]},
With[{notAppears=Table[Union[Flatten[appears]],{Length[eqns]}]},
With[{vPairs=ReplaceAll[appears,xx_[tt_]:>{First[First[Position[sv,xx]]],tt-Global`t+1}/;xx=!=List],
theSubs=Thread[allv->varsForSubbing]},
With[{listOfFuncs=Map[Apply[Function,#]&,
MapThread[{#1/.theSubs,({#2}/.theSubs)}&,{appears,eqns}]]},
With[{rawDrvLists=Map[Flatten,
MapThread[doDrvsToSecond,{listOfFuncs,appears}]]},
With[{drvLists=rawDrvLists//.linPtSubs,
dListOfFuncs=Map[Apply[Function,#]&,
MapThread[{#1/.theSubs,({#2}/.theSubs)}&,{appears,rawDrvLists}]]},
If[linPtSubs==={},{{},sv,appears,listOfFuncs,sv,vPairs,drvLists},
{linPtSubs,sv,appears,listOfFuncs,sv,vPairs,
drvLists,rawDrvLists,dListOfFuncs}]]]]]]]]]])/;And[
AMAEquations[model]=!=AMAundefinedModel,
FreeQ[hmatSSSolnSubs[model,opts],FindRoot]]
makeMatlabFuncs[___]:=AMAundefinedModel;


constructFDerivs[theModel_Symbol,opts___?OptionQ]:=
With[{mFuncStuff=
makeMatlabFuncs[theModel,opts],
cut=
errorVarCutoff[theModel],
neq=AMANeq[theModel],
nlags=AMALags[theModel]},
With[{theFuncs=mFuncStuff[[4]],
theVars=mFuncStuff[[3]],
ssSubs=mFuncStuff[[1]],
ssVals=
AMASS[theModel,opts],
vAppByFunc=mFuncStuff[[6]],
vApp=
variablesAppearing[mFuncStuff[[6]],cut]},
With[{theFuncsApplied=
Flatten[MapThread[Apply[#1,#2]&,{theFuncs,theVars}]],
vMappers=Map[variableMapping[vApp,#]&,vAppByFunc],
modNz=constructNz[vApp,neq,nlags],
modLeads=constructLeads[vApp],
modEta=constructEta[vApp,cut]},
With[{firstDrvs=
(MapThread[doFacet[1,#1,#2]&,{theFuncsApplied,theVars}]/.ssSubs)/.
ssVals},
With[{miVals=Map[miInFacet[Length[#],2]&,theVars]},
With[{secondDrvs=
(MapThread[doFacetMi,
{theFuncs,theVars,miVals}]/.ssSubs)/.ssVals},
With[{jaggedFuncVals=
Map[Apply[Join,#]&,
Transpose[{firstDrvs,secondDrvs}]],
squLen=MultiIndex`pascal[Length[vApp],2]-1,
forExplode=Map[
Function[xx,Apply[variableExplode[#1,Length[#2]]&,xx]],vMappers]},
{doExplosions[jaggedFuncVals,forExplode,squLen],
modNz,modLeads,modEta,nlags}
]]]]]]]



doAnExplosion[lilDerivs_List,targets_List,squareDim_Integer]:=
Module[{resMat=Table[0,{squareDim}]},
resMat[[targets]]=lilDerivs;
resMat]

doExplosions[jaggedDerivs_List,targets_List,squareDim_Integer]:=
MapThread[doAnExplosion[#1,#2,squareDim]&,
{jaggedDerivs,targets}]


variableMapping[longVec_List,shortVec_List]:=
With[{varNos=Flatten[Map[Position[longVec,#]&,shortVec]]},
With[{mapper=Function[xx,
Module[{resVec=Table[0,{Length[longVec]}]},
resVec[[varNos]]=xx;resVec]]},
{mapper,varNos}]]

numInFacet[dim_Integer,deg_Integer]:=
(MultiIndex`pascal[dim,deg]-MultiIndex`pascal[dim,deg-1])/;And[dim>0,deg>0]
numInFacet[dim_Integer,0]:=1/;dim>0

doFacetMi[aFunc_,theVars_,theMi_List]:=Flatten[
Map[Apply[(Apply[Derivative,#][aFunc]),theVars]&,theMi]]


miInFacet[dim_Integer,deg_Integer]:=
Drop[NestList[MultiIndex`algorithm2,
predMI[dim,deg],numInFacet[dim,deg]],1]/;And[dim>0,deg>0]



theIndex[beta_List]:=index[beta];(*for export*)


index[beta_List]:=
With[{nn=Length[beta]},
Sum[With[{sc=(Apply[Plus , Drop[beta,cc-1]])-1},
MultiIndex`pascal[nn-cc+1,sc]],{cc,nn}]]




Unprotect[Power]
0^0=1
0.^0=1
Protect[Power]



pyramidPowersJava[pyramidVals_?MatrixQ,ord_Integer]:=
Module[{},
If[MatrixQ[pyramidVals,NumberQ],
Module[{},
With[{numPyrs=Length[pyramidVals],pyrCols=Length[pyramidVals[[1]]]},
timeResult["javapyramidpowerscall:",
Partition[MultiIndex`pyramidPowers[
Flatten[pyramidVals],ord,numPyrs],
pyrCols]]]],
pyramidPowers[pyramidVals,ord]]]




pyramidPowers[pyramidVals_?MatrixQ,ord_Integer]:=
With[{numPyrs=Length[pyramidVals]},
With[{preLam=Table[0,{numPyrs}]},
With[{lamList=Drop[NestList[MultiIndex`algorithm2,preLam,
MultiIndex`pascal[ord,numPyrs]-1],1]},If[MatrixQ[pyramidVals,NumberQ],
(*sparse array does't seem to help much for numbers here*)
SparseArray[cOuter[lamList,Transpose[pyramidVals]]],
(*sparse array very slow for symbolic here*)
Outer[powerFunc,lamList,Transpose[pyramidVals],1,1]]]]]




cOuter=Compile[{{lamList,_Integer,2},{pyramidVals,_Real,2}},
Module[{},
With[{interAct=Flatten[Outer[List,lamList,pyramidVals,1,1],1]},
Partition[Flatten[
Map[appCPowerFunc,interAct]],Length[pyramidVals]]]]]


appCPowerFunc[{lamVals_?VectorQ,pyramidVals_?VectorQ}]:=
Module[{},
cPowerFunc[Round[lamVals],pyramidVals]]




powerFunc[kVals_?VectorQ,pyrVals_?VectorQ]:=
Module[{},
With[{prs=DeleteCases[Transpose[{pyrVals,kVals}],{_,0}]},
If[Or[prs==={},MemberQ[prs,{0.0,_}|{0,_}]],0,
With[{unPrs=Transpose[prs]},
If[VectorQ[unPrs[[1]],NumberQ],
cPowerFunc[unPrs[[2]],unPrs[[1]]][[1,1]],
Apply[Times,unPrs[[1]]^unPrs[[2]]]]]]]]


cPowerFunc=
Compile[{{kVals,_Integer,1},{pyrVals,_Real,1}},
Module[{},
With[{prs=Transpose[{pyrVals,kVals}]},
With[{selZap=Apply[Plus,Map[If[And[#[[1]]==0,#[[2]]!=0],1,0]&,prs]]},
If[selZap>0,{{0}},
With[{relevant=Transpose[Select[prs,And[#[[1]]!=0,#[[2]]!=0]&]]},
With[{res=Apply[Times,(relevant[[1]]^Round[relevant[[2]]])]},
{{res}}]]]]]],{{res,_Real},{Count[___],_Integer}}]



doFacet[deg_Integer,aFunc_,theVars_List]:=
If[deg==1,
Map[D[aFunc,#]&,theVars]]

predMI[dim_Integer,deg_Integer]:=
Join[Table[0,{dim-1}],{deg-1}]/;And[dim>0,deg>0]

gMults02[gPyr01_?MatrixQ,
assumedMomSubs_List:{}]:=
Module[{resultMat},
With[{deg=2,
ss=numCols[gPyr01],
bigNu=numRows[gPyr01],
pPows=
pyramidPowersJava[gPyr01,2]},
With[{preNu=Table[0,{numCols[gPyr01]}],
nuDimInfo=dim02Cols[ss],
lamDimInfo=dim02Cols[bigNu]},
With[{nuTwo=nuDimInfo[[1]],nuNotTwo=nuDimInfo[[2]],
lamTwo=lamDimInfo[[1]],lamNotTwo=lamDimInfo[[2]]},
With[{
gm22=
(pPows[[lamTwo+bigNu,All]]),
gm21=
(2*pPows[[lamNotTwo+bigNu,All]])},
resultMat=
zeroMatrix[numInFacet[bigNu,2],numInFacet[ss,2]];
If[ss>1,
With[{
gm12=
Transpose[Map[(pPows[[Range[bigNu],#[[1]]]]*
pPows[[Range[bigNu],#[[2]]]])&,
nuDim1LamDim2Case[ss]]],
gm11=
Map[compCaseFunc[pPows],
nuDim1LamDim1Case[ss,bigNu],{2}]
},
resultMat[[lamTwo,nuTwo]]=gm22;
resultMat[[lamNotTwo,nuTwo]]=gm21;
resultMat[[lamTwo,nuNotTwo]]=gm12;
resultMat[[lamNotTwo,nuNotTwo]]=gm11];,
resultMat[[lamTwo,nuTwo]]=gm22;
resultMat[[lamNotTwo,nuTwo]]=gm21];
resultMat
]]]]]


compCaseFunc[xx_List]:=(Function[yy,(xx[[yy[[1]],yy[[3]]]]*
xx[[yy[[2]],yy[[4]]]]+
xx[[yy[[2]],yy[[3]]]]*
xx[[yy[[1]],yy[[4]]]])])


nuDim1LamDim2Case[miLen_Integer]:=
With[{chunks=Range[miLen-1,1,-1]},Transpose[{
Flatten[Table[Table[ii,{jj,chunks[[ii]]}],{ii,miLen-1}]],
Flatten[Table[Table[ii+jj,{jj,chunks[[ii]]}],{ii,miLen-1}]]}]
]




nuDim1LamDim1Case[miLen_Integer,bigNu_Integer]:=
With[{n1l2=nuDim1LamDim2Case[miLen],
forLam=nuDim1LamDim2Case[bigNu]},
Flatten[Outer[Join,forLam,n1l2,1],0]]



dim02Cols[miLen_Integer]:=
If[miLen==1,{{1},{}},
With[{gaps=Range[miLen,2,-1],
allVals=numInFacet[miLen,2]},
With[{yesVals=FoldList[Plus,1,gaps]},
{yesVals,Complement[Range[allVals],yesVals]}]]]


gMults02Java[gPyr01_?MatrixQ,
assumedMomSubs_List:{}]:=
If[MatrixQ[gPyr01,NumberQ],
Module[{},
With[{numPyrs=numRows[gPyr01],pyrCols=numCols[gPyr01]},
Partition[MultiIndex`gMults02[
Flatten[gPyr01],numPyrs],
numInFacet[pyrCols,2]]]],gMults02[gPyr01,assumedMomSubs]]


gMults02Java[gPyr01_?MatrixQ,
assumedMomSubs_List:{}]:=
If[MatrixQ[gPyr01,NumberQ],
Module[{},
With[{numPyrs=numRows[gPyr01],pyrCols=numCols[gPyr01]},
Partition[MultiIndex`gMults02[
Flatten[gPyr01],numPyrs],
numInFacet[pyrCols,2]]]],gMults02[gPyr01,assumedMomSubs]]


(*<<fromExperPattern.mth*)
gMults03[gPyr02_?MatrixQ,gArgs_Integer,
assumedMomSubs_List:{}]:=MulMuGMult[02,gArgs,gPyr02][[1]]


notGMults02[gPyr01_?MatrixQ,
assumedMomSubs_List:{}]:=
Module[{resultMat},
With[{deg=2,
ss=numCols[gPyr01],
bigNu=numRows[gPyr01],
pPows=
pyramidPowersJava[gPyr01,2]},
With[{preNu=Table[0,{numCols[gPyr01]}],
nuDimInfo=dim02Cols[ss],
lamDimInfo=dim02Cols[bigNu]},
With[{nuTwo=nuDimInfo[[1]],nuNotTwo=nuDimInfo[[2]],
lamTwo=lamDimInfo[[1]],lamNotTwo=lamDimInfo[[2]]},
With[{
gm22=
(pPows[[lamTwo+bigNu,All]]),
gm21=
(2*pPows[[lamNotTwo+bigNu,All]])},
resultMat=
zeroMatrix[numInFacet[bigNu,2],numInFacet[ss,2]];
If[ss>1,
With[{
gm12=
Transpose[Map[(pPows[[Range[bigNu],#[[1]]]]*
pPows[[Range[bigNu],#[[2]]]])&,
nuDim1LamDim2Case[ss]]],
gm11=
Map[compCaseFunc[pPows],
nuDim1LamDim1Case[ss,bigNu],{2}]
},
resultMat[[lamTwo,nuTwo]]=gm22;
resultMat[[lamNotTwo,nuTwo]]=gm21;
resultMat[[lamTwo,nuNotTwo]]=gm12;
resultMat[[lamNotTwo,nuNotTwo]]=gm11];,
resultMat[[lamTwo,nuTwo]]=gm22;
resultMat[[lamNotTwo,nuTwo]]=gm21];
resultMat
]]]]]


AMASS[model_,opts___?OptionQ]:=
(If[FreeQ[Flatten[{opts}],hmatLinPtSubs],
computeSteadyState[model,opts],With[{
hres=ReplaceAll[hmatLinPtSubs,Flatten[{opts}]]},
If[hres==={},computeSteadyState[model,opts],hres]]])/;
AMAEquations[model]=!=AMAundefinedModel;
AMASS[___]:=AMAundefinedModel;


AMASSGuess[_Symbol]:={initSSGuess->{}}
Options[computeSteadyState]={initSSGuess->{}};


computeSteadyState[model_Symbol,opts___?OptionQ]:=
model/:computeSteadyState[model,opts]=
With[{eqns=AMAEquations[model]},
With[{epsVars=Union[Cases[eqns,Global`eps[x_][_]->Global`eps[x],Infinity]]},
With[{sv=Complement[stateVariables[eqns],epsVars]},
With[{initGuess=ReplaceAll[(ReplaceAll[ReplaceAll[ReplaceAll[
initSSGuess,Flatten[{opts}]],AMASSGuess[model]],
   {initSSGuess->Table[0,{Length[sv]}]}]),{{}->Table[0,{Length[sv]}]}],
exogXform=ReplaceAll[ReplaceAll[exogenizeVars,Flatten[{opts}]],
{exogenizeVars->{}}]
},
With[{ssVars=Transpose[{ReplaceAll[
Through[sv[Global`t]],ssSubs[]],initGuess}],
ssEqns=Thread[(ReplaceAll[eqns,ssSubs[]])==0]},
With[{theSoln=
Apply[FindRoot , (ReplaceAll[Join[{ssEqns},ssVars],exogXform])]},
If[FreeQ[theSoln,FindRoot],
newAMAModel[model,AMAEquations[model],(initSSGuess->Map[Last,theSoln])]];
theSoln]]]]]]


ssSubs[]:=
{Global`eps[_][_]->0,x_Symbol[Global`t+_.]:>Symbol[SymbolName[x]<>"AMAss"]};




hmatSSSolnSubs[model_Symbol,opts___?OptionQ]:=
(model/:ssSolnSubs[model,opts]=
{ssSubs[],AMASS[model,opts]})/;AMAEquations[model]=!=AMAundefinedModel;
hmatSSSolnSubs[___]:=AMAundefinedModel;



genericHmatSSSolnSubs[model_Symbol,opts___?OptionQ]:=
(model/:ssSolnSubs[model,opts]=
{ssSubs[]})/;AMAEquations[model]=!=AMAundefinedModel;
genericHmatSSSolnSubs[___]:=AMAundefinedModel;


conStochBs[theModel_Symbol,augmentedSubs_List:{}]:=
With[{theNewB=getNewB[theModel],
theBMat=AMABmat[theModel],
theSMat=-AMAS0Inv[theModel] . AMATheta[theModel],
epsVars=getEpsVars[theModel],
numErrs=getNumErrs[theModel],
leads=getLeadsNeeded[theModel],
etaM=getEtaM[theModel],
theNz=getNzCols[theModel],neq=getEta0[theModel],nlags=getNLags[theModel]},
With[{maxLead=Max[Map[Last,getLeadsNeeded[theModel]]],
neqRange=etaM+Range[neq],
tkVals=toKeep[theNz,neq,nlags,numErrs]},
With[{sigCol=theNewB[[neqRange,etaM+{1}]]},
With[{iterBMat=makeIterBMat[
neq,nlags,numErrs,theBMat,theSMat,sigCol[[Range[neq]]]]
},
With[{},
With[{bigMats=
Map[First,NestList[{stochAugB[iterBMat . #[[1]],epsVars,#[[2]]],#[[2]]+1}&,
{IdentityMatrix[numCols[iterBMat]],1},maxLead+1]]},
With[{keeperMats=Map[#[[tkVals,tkVals]]&,bigMats]},
getStochConBs[theModel]^=keeperMats;
keeperMats
]]]]]]]


toKeep[theNz_List,neq_Integer,nlags_Integer,numErrs_Integer]:=
With[{lagMoreOne=Range[1,neq*(nlags-1)]},
Union[Intersection[theNz,lagMoreOne],neq*(nlags-1)+Range[neq+numErrs+1]]]


laggedToKeep[theNz_List,neq_Integer,nlags_Integer]:=
Range[Length[Select[theNz,#<=neq*(nlags-1)&]]]


nonLaggedToKeep[theNz_List,neq_Integer,nlags_Integer]:=
(Select[theNz,#>neq*(nlags-1)&]-neq*(nlags-1))

makeIterBMat[neq_Integer,nlags_Integer,numErrs_Integer,
theBMat_?MatrixQ,theSMat_?MatrixQ,theSigMat_?MatrixQ]:=
With[{bPart=If[nlags==1,theBMat,
blockMatrix[{{Drop[IdentityMatrix[neq*nlags],neq]},{theBMat}}]],
sPart=If[nlags==1,theSMat,
blockMatrix[{{zeroMatrix[neq*(nlags-1),numErrs]},{theSMat}}]],
sigPart=If[nlags==1,theSigMat,
blockMatrix[{{zeroMatrix[neq*(nlags-1),1]},{theSigMat}}]]
},
blockMatrix[{{bPart,sPart,sigPart}}]]


makeSmallIterBMat[neq_Integer,nlags_Integer,numErrs_Integer,
theBMat_?MatrixQ,theSMat_?MatrixQ,theSigMat_?MatrixQ]:=
With[{bPart=If[nlags==1,theBMat,
blockMatrix[{{Drop[IdentityMatrix[neq*nlags],neq]},{theBMat}}]],
sPart=If[nlags==1,theSMat,
blockMatrix[{{zeroMatrix[neq*(nlags-1),numErrs]},{theSMat}}]],
sigPart=If[nlags==1,theSigMat,
blockMatrix[{{zeroMatrix[neq*(nlags-1),1]},{theSigMat}}]]
},
numericEliminateInessentialLags[
If[nlags==1,
{blockMatrix[{{bPart,sPart,sigPart}}],Range[neq*nlags]},
{blockMatrix[{{blockMatrix[{{zeroMatrix[neq*(nlags-1),neq],IdentityMatrix[neq*(nlags-1)]}}]},(*probsbly*)
{blockMatrix[{{bPart,sPart,sigPart}}]}}],Range[neq*nlags]}]]]


getSBDerivs[___]={}



eqStochSysCon01[theModel_?isAMAModel,assumedMomSubs_List:{}]:=
With[{sBDerivs=getSBDerivs[theModel,assumedMomSubs]},
If[Length[sBDerivs]>0,sBDerivs[[1]],
With[{
bDerivs=getBDerivs[theModel],
modNz=getNzCols[theModel],
leadsNeeded=getLeadsNeeded[theModel],
nlags=getNLags[theModel],
modEta={getEtaM[theModel],eta0=getEta0[theModel],etaP=getEtaP[theModel]},
sFDerivs=getFDerivs02[theModel],
numErrs=getNumErrs[theModel],
augmentedSubs=Join[makeSomeMomSubs[getEpsVars[theModel]],assumedMomSubs]
},
With[{foStochRes=doFirstOrderStoch[sFDerivs,bDerivs,modEta,
modNz,nlags,getEpsVars[theModel],leadsNeeded,
augmentedSubs]},
theModel /: getSBDerivs[theModel,assumedMomSubs]=
{Drop[foStochRes[[4]],getEtaM[theModel]]};
getNewB[theModel]^=foStochRes[[4]];
foStochRes
]]]]




doFirstOrderStoch[sFDerivs_List,bDerivs_?MatrixQ,
theEta:{etaM_Integer,neq_Integer,etaP_Integer},modNz_List,nlags_Integer,
epsVars_List,leadsNeeded_List,
assumedMomSubs_List:{}]:=
Module[{theBMat=zeroMatrix[neq,neq*nlags],
numErrs=Length[epsVars]},
theBMat[[All,modNz]]=bDerivs[[etaM+Range[neq],Range[Length[modNz]]]];
With[{
theS0Mat=bDerivs[[etaM+Range[neq],Length[modNz]+Range[numErrs]]]},
With[{sigN001=cmpSig01[sFDerivs,
{etaM,neq,etaP},modNz,theBMat,theS0Mat,
epsVars,leadsNeeded,assumedMomSubs]},
{theBMat,theS0Mat,sigN001,
blockMatrix[{{bDerivs[[All,Range[etaM]]],
blockMatrix[{{zeroMatrix[etaM,1]},{sigN001}}]}}]}]]]


cmpSig01[sFDerivs_List,theEta:{etaM_Integer,eta0_Integer,etaP_Integer},
theNz_List,
bMat_?MatrixQ,s0Mat_?MatrixQ,epsVars_List,theLeads_List,
assumedMomSubs_List:{}]:=
With[{(*eta0=getEta0[theEta],*)
nleads=Max[Map[Last,theLeads]],
someMomSubs=Join[makeSomeMomSubs[epsVars],assumedMomSubs]},
With[{nlags=numCols[bMat]/eta0},
With[{extendedBMat=
multExtndMult[bMat,initPath[bMat,nlags],nleads-1]},
With[{
initMat=If[nlags==1,IdentityMatrix[eta0],
blockMatrix[{{zeroMatrix[eta0,(nlags-1)*eta0],
IdentityMatrix[eta0]}}]]},
With[{sigMults=blockMatrix[Transpose[(*probably*)
{FoldList[Plus,
initMat,Drop[Partition[Drop[extendedBMat,eta0*(nlags-1)],eta0],0]]}]]},
With[{forMult=conformToF0P[sigMults,theEta,theLeads],
numErrs=etaM-Length[theNz]},
(*sparse arrays don't seem to help much for sedmodelopt*)
With[{lud=LinearSolve[Identity[
sFDerivs[[1]][[All,etaM+Range[eta0+etaP]]]].
Identity[Chop[Expand[forMult]]]]},
With[{epsTimeList=conformToFP[
fleshOutErrs[epsVars,s0Mat,bMat,nleads],theEta,theLeads]},
(*need this expand sub or error*)
With[{},With[{ludApplied=Chop[Expand[lud[Normal[
(Chop[Expand[-sFDerivs[[1]][[All,etaM+eta0+Range[etaP]]] . 
epsTimeList]//.someMomSubs])]]]]},
forMult. ludApplied]]]]]]]]]]

fleshOutErrs[epsVars_List,s0Mat_?MatrixQ,bMat_?MatrixQ,nLeads_Integer]:=
With[{neq=numRows[bMat],numErrs=Length[epsVars]},
With[{nlags=numCols[bMat]/neq},
With[{bigS0Mat=If[nlags==1,
s0Mat,blockMatrix[{{
zeroMatrix[neq*(nlags-1),numErrs]},{s0Mat}}]]},
With[{impact=Map[bigS0Mat.#&,
Map[Transpose[{Through[epsVars[t+#]]}]&,
Range[nLeads]]]},
With[{theMap=MapIndexed[Join[
zeroMatrix[neq*(#2[[1]]-1),
1],
multExtndMult[bMat,#1,nLeads-#2[[1]]]]&,
impact]},
With[{doBMult=Apply[Plus,
theMap]},
doBMult]]]]]]



initPath[origB_?MatrixQ,nlags_Integer]:=
With[{eta0=numRows[origB]},
If[nlags==1,origB,
blockMatrix[{{Drop[IdentityMatrix[eta0*nlags],eta0]},{origB}}]]]



conformToF0P[ceDrvs01_?MatrixQ,
{etaM_Integer,eta0_Integer,etaP_Integer},leadsNeeded_List]:=
With[{leadPos=Flatten[Join[Map[leadToPos[#,eta0]&,leadsNeeded]]],
rightCols=Range[eta0]+numCols[ceDrvs01]-eta0},
ceDrvs01[[Join[Range[eta0],eta0+leadPos],rightCols]]]


conformToFP[ceDrvs01_?MatrixQ,
{etaM_Integer,eta0_Integer,etaP_Integer},leadsNeeded_List]:=
With[{leadPos=Flatten[Join[Map[leadToPos[#,eta0]&,leadsNeeded]]]},
ceDrvs01[[leadPos]]]


leadToPos[{varNo_Integer,leadVal_Integer},neq_Integer]:=(leadVal-1)*neq+varNo




AMAHmat[model_,opts___?OptionQ]:=
(model/:AMAHmat[model,opts]=makeHmat[model,opts])/;
AMAEquations[model]=!=AMAundefinedModel;
AMAHmat[___]:=AMAundefinedModel;



makeHmat[model_Symbol,opts___?OptionQ]:=
(With[{eqns=AMAEquations[model]},
With[{linPtSubs=Flatten[hmatSSSolnSubs[model,opts]]},
With[{allv=DeleteCases[completeArgList[eqns],Global`eps[_][_]]},
With[{raw=Outer[D,eqns,allv]},
If[linPtSubs==={},{{},raw},{ReplaceRepeated[raw,linPtSubs],raw}]]]]])/;And[
AMAEquations[model]=!=AMAundefinedModel,
FreeQ[hmatSSSolnSubs[model,opts],FindRoot]]

makeHmat[model_Symbol,opts___?OptionQ]:=
(With[{eqns=AMAEquations[model]},
With[{linPtSubs=Flatten[genericHmatSSSolnSubs[model,opts]]},
With[{allv=DeleteCases[completeArgList[eqns],Global`eps[_][_]]},
With[{raw=Outer[D,eqns,allv]},
If[linPtSubs==={},{{},raw},{ReplaceRepeated[raw,linPtSubs],raw}]]]]])/;And[
AMAEquations[model]=!=AMAundefinedModel]
makeHmat[___]:=AMAundefinedModel;







AMASoln[model_,opts___?OptionQ]:=
With[{finalResultForTiming=
(model/:AMASoln[model,opts]=
With[{ilag=ReplaceAll[ReplaceAll[infoLag,Flatten[{opts}]],infoLag->0]},
Chop[numericAMA[AMAHmat[model,opts][[1]],AMALags[
model],ilag]]])},
finalResultForTiming]/;VectorQ[Flatten[AMAHmat[model,opts][[1]]],NumberQ]





AMABmat[model_,opts___?OptionQ]:=(AMASoln[model,opts][[-3]])/;
AMAEquations[model]=!=AMAundefinedModel;
AMABmat[___]:=AMAundefinedModel;
AMAQmat[model_,opts___?OptionQ]:=(AMASoln[model,opts][[-4]])/;
AMAEquations[model]=!=AMAundefinedModel;
AMAQmat[___]:=AMAundefinedModel;
AMASmat[model_,opts___?OptionQ]:=(AMASoln[model,opts][[-2]])/;
AMAEquations[model]=!=AMAundefinedModel;
AMASmat[___]:=AMAundefinedModel;
AMAS0Inv[model_,opts___?OptionQ]:=(AMASoln[model,opts][[-1]])/;
AMAEquations[model]=!=AMAundefinedModel;
AMAS0Inv[___]:=AMAundefinedModel;
AMAEigSystem[model_,opts___?OptionQ]:=(AMASoln[model,opts][[-5]])/;
AMAEquations[model]=!=AMAundefinedModel;
AMAEigSystem[___]:=AMAundefinedModel;

AMATheta[model_,opts___?OptionQ]:=
(model/:AMATheta[model,opts]=
With[{eqns=AMAEquations[model]},
With[{epsVars=Union[Cases[eqns,Global`eps[x_][_]->Global`eps[x],Infinity]]},
With[{theTheta=Outer[D,eqns,Through[epsVars[Global`t]]]},
theTheta]]])



cmpCESolution02[
{fMat_?MatrixQ,y0InvPart_,interMedPart_?MatrixQ},theConstant_?VectorQ,
theModel_Symbol]:=
With[{
eta0=getEta0[theModel],
etaP=getEtaP[theModel]},
With[{backwardLookingSoln=
(sFyDoKronLUBack[eta0,y0InvPart,theConstant])},
With[{secTerm=
SparseArray[backwardLookingSoln]},
With[{firTerm=
SparseArray[(interMedPart)]},
With[{interMed=
firTerm. secTerm},
With[{bFuncDrvSoln=
backwardLookingSoln-interMed},
With[{matSoln=
Transpose[Partition[
Flatten[bFuncDrvSoln],eta0]],
theFuture=
Transpose[Partition[Flatten[fMat. bFuncDrvSoln],etaP]]},
With[{fullMatSoln=
blockMatrix[{{matSoln},{theFuture}}]},
fullMatSoln
]]]]]]]]






cmpStochWithCE02[theModel_Symbol,ceSoln_?MatrixQ,ceSels_List,stochSels_List,
{stochF0Inv_LinearSolveFunction,stochMidBeg_?MatrixQ},theConstant_?VectorQ,theLinSys_List,
assumedMomSubs_List:{}]:=
With[{etaP=getEtaP[theModel],eta0=getEta0[theModel],etaM=getEtaM[theModel]},
With[{yp=(getFDerivs02[theModel][[1]])[[All,Range[etaP]+etaM+eta0]]},
With[{
ceDiff=kronI[numInFacet[getEtaM[theModel]+1,2],yp][[stochSels]] .
theLinSys[[-1]][[All,ceSels]] . Flatten[Transpose[ceSoln[[Range[eta0]]]]]},
With[{theNoFuture=sFyDoKronLUBack[getEta0[theModel],
stochF0Inv,
Flatten[theConstant[[stochSels]]-ceDiff]]},
Transpose[Partition[Flatten[
(theNoFuture -stochMidBeg . theNoFuture)],eta0]]/.assumedMomSubs]]]]





sFyDoKronLUBack[neq_Integer,luD_LinearSolveFunction,rhs_?VectorQ]:=
Module[{resMat=zeroMatrix[Length[rhs],1]},
With[{otherDim=Length[rhs]/neq},
Do[
With[{theRes=(
luD[rhs[[(ii-1)*neq+Range[neq]]]])},
(resMat[[Range[neq]+(ii-1)*neq,1]]=theRes)],
{ii,otherDim}];
resMat]]



cmpStochConstant02[theModel_Symbol,assumedMomSubs_List:{}]:=
With[{
augmentedSubs=Join[makeSomeMomSubs[getEpsVars[theModel]],assumedMomSubs],
leadsNeeded=getLeadsNeeded[theModel],
modNz=getNzCols[theModel],
sFDerivs=getFDerivs02[theModel],
theStochBs=getStochConBs[theModel],
etaM=getEtaM[theModel],
eta0=getEta0[theModel],
etaP=getEtaP[theModel],
nlags=getNLags[theModel]
},
With[{
laggedVarsRequired=laggedToKeep[modNz,eta0,nlags],
stateVarDrvs=blockMatrix[{{IdentityMatrix[etaM],zeroMatrix[etaM,1]}}],
bFuncsRequired=stochForMult[leadsNeeded,eta0],
bColsRequired=stochSelectMin[{etaM,eta0,etaP},modNz,nlags]
},
With[{
bAndBCompDrvs=
Flatten[Drop[theStochBs,1][[All,
Length[laggedVarsRequired]+Range[eta0],bColsRequired]],1]},
With[{gFuncDrvs=blockMatrix[{{stateVarDrvs},
{bAndBCompDrvs[[bFuncsRequired,All]]}}]},
(*needs subs or else big time hit for sedmodelopt*)
With[{theGMultVals=Expand[gMults02[gFuncDrvs,augmentedSubs]]//.augmentedSubs},
With[{theConstant=(-1)*sFDerivs[[2]] . theGMultVals},
Flatten[Transpose[theConstant]]
]]]]]]


doBPows[stochBs_List,{etaM_Integer,neq_Integer,etaP_Integer},
modNz_List,nlags_Integer]:=
With[{ltk= laggedToKeep[modNz,eta0,nlags]},
Print["dims",{Dimensions[stochBs],{etaM,neq,etaP}}];
Map[nullToNullSparseArray[#[[Length[ltk]+Range[neq],Length[ltk]+Range[neq]]]]&,stochBs]]



stochKronSum02[stochBs_List,{etaM_Integer,neq_Integer,etaP_Integer},
modNz_List,nlags_Integer,lilRel_List,assumedMomSubs_List:{}]:=
With[{supRows=relRowsSplitUp02[lilRel,neq,Length[stochBs]*neq]},
With[{ssm= stochSelectMin[{etaM,neq,etaP},modNz,nlags]},
With[{bPowsLil=doBPows[stochBs,{etaM,neq,etaP},modNz,nlags],
lilStochBs=Map[#[[ssm,ssm]]&,stochBs]},
(* need expand sub or error*)
With[{multListLil=
Map[SparseArray[Expand[gMults02Java[#,assumedMomSubs]]//.assumedMomSubs]&,
Drop[lilStochBs,1]]},
With[{theINuLil=spIdentityMatrix[numRows[multListLil[[1]]]]},
With[{},
With[{
specRes=SparseArray[Flatten[
MapIndexed[specificKronSum[#1,#2[[1]],
multListLil,bPowsLil]&,supRows],1]]},
{specRes}]]]]]]]

nullToNullSparseArray[aMat_?MatrixQ]:=
If[aMat==={},{},SparseArray[aMat]]




cmpStochLinSys02[theModel_Symbol,augmentedSubs_List]:=
With[{
sFDerivs=getFDerivs02[theModel],
theStochBs=getStochConBs[theModel],
modEta=getEta[theModel],
nlags=getNLags[theModel],
neq=getEta0[theModel],
leadsNeeded=getLeadsNeeded[theModel],
modNz=getNzCols[theModel],
stateVarDim=getEtaM[theModel]+1},
With[{drvsInFacet=numInFacet[stateVarDim,2]},
With[{(*relevantRows=stochKeepers[leadsNeeded,neq,drvsInFacet],*)
lilRel=stochKeepers[leadsNeeded,neq,1]},
With[{kronSums=
 stochKronSum02[Drop[Expand[theStochBs],
-1],modEta,modNz,nlags,lilRel,augmentedSubs]},
With[{fMat=kronSums[[1]][[
reOrderCols[leadsNeeded,drvsInFacet]]]},
Join[{fMat}]]]]]]





specificKronSum[spRows_List,theLead_Integer,multVals_List,bPows_List]:=
If[spRows==={},{},
With[{nuDim=Length[multVals[[1]]]},
With[{theINuLil=spIdentityMatrix[nuDim]},
With[{anchor=kron[theINuLil,bPows[[theLead+1]][[spRows]]]},
With[{guts=Map[Apply[kron,rightBWPair[spRows,bPows[[Range[theLead]]],multVals,#1]
(*Transpose[SparseArray[multVals[[#1]]]],
SparseArray[bPows[[Range[theLead]]][[-#1]]][[spRows]]
*)]&,
Range[theLead]]},
With[{theIter=(Apply[Plus,guts])},
anchor+theIter
]]]]]]

rightBWPair[spRows_List,bVals_List,wVals_List,theIndx_Integer]:=
{Transpose[SparseArray[wVals[[theIndx]]]],
SparseArray[bVals][[-theIndx]][[spRows]]}




reOrderCols[leadsNeeded_List,nuDim_Integer]:=
With[{maxLead=Max[Map[Last,leadsNeeded]]},
With[{regrouped=Map[ofGivenLead[leadsNeeded,#]&,Range[maxLead]]},
With[{theEtaKs=Map[Length,regrouped]},
With[{theSums=Drop[FoldList[Plus,0,theEtaKs],-1]},
With[{forEachDrv=DeleteCases[MapThread[rangeForK,{theEtaKs,nuDim*theSums}],{}]},
With[{allOfEm=Flatten[
Map[forEachDrv+bumpByEtaK[DeleteCases[theEtaKs,0],#]&,Range[nuDim]]]},
allOfEm]]]]]]



stochForMult[leadsNeeded_List,neq_Integer]:=
Union[Range[neq],Map[aLeadPos[#,neq]&,leadsNeeded]]
aLeadPos[{varNo_Integer,aLead_Integer},neq_Integer]:=
neq+(aLead-1)*neq+varNo


stochKeepers[leadsNeeded_List,stochNeq_Integer,drvDim_Integer]:=
With[{rightOrd=Sort[leadsNeeded,cmpByLeadVarNo]},
Flatten[Transpose[
Map[oneStochKeeper[#[[1]],#[[2]],stochNeq,drvDim]&,leadsNeeded]]]]


cmpByLeadVarNo[{varNoA_Integer,theLeadA_Integer},
{varNoB_Integer,theLeadB_Integer}]:=
Or[theLeadA<theLeadB,And[theLeadA==theLeadB,varNoA<varNoB]]



oneStochKeeper[varNo_Integer,theLead_Integer,stochNeq_Integer,drvDim_Integer]:=
With[{varPos=varNo+Range[0,drvDim*stochNeq-1,stochNeq]},
varPos+(theLead-1)*stochNeq*drvDim]




stochSelectMin[modEta:{etaM_Integer,eta0_Integer,etaP_Integer},
modNz_List,nlags_Integer]:=
With[{(*eta0=getEta0[modEta],etaM=getEtaM[modEta],*)numNonErrs=Length[modNz]},
With[{numErrs=etaM-numNonErrs,
ltk=laggedToKeep[modNz,eta0,nlags],
nltk=nonLaggedToKeep[modNz,eta0,nlags]},
With[{lenLtk=Length[ltk]},
Join[Range[lenLtk],lenLtk+nltk,lenLtk+eta0+Range[1,numErrs+1]]]]]




sFyCmpY0YPInv[sFDerivs_List,theEta:{etaM_Integer,eta0_Integer,etaP_Integer},
fVars_Integer,fMat_?MatrixQ]:=
With[{numDrvs=numCols[fMat]/fVars},
With[{f0Raw=sFDerivs[[1]][[All,etaM+Range[eta0]]],
fpRaw=sFDerivs[[1]][[All,etaM+eta0+Range[etaP]]]},
With[{f0Inv=
LinearSolve[Normal[f0Raw]]},
With[{f0IFP=
f0Inv[fpRaw]},
With[{kF0IFP=
kronI[numDrvs,f0IFP]},
With[{midTerm=
LinearSolve[Normal[
spIdentityMatrix[numDrvs*etaP]+Expand[SparseArray[
fMat .kF0IFP]]]]},
With[{curious=Map[(
SparseArray[Chop[kF0IFP.
midTerm[Normal[#]]]])&,Transpose[fMat]]},
With[{midBeg=
kF0IFP.
midTerm[fMat]},
{f0Inv,midBeg}]]]]]]]]

sFyCmpY0YP[fDerivs_List,{etaM_Integer,eta0_Integer,etaP_Integer},
numVars_Integer]:=
With[{numDrvs=numInFacet[numVars,2]},
With[{f0Raw=fDerivs[[All,etaM+Range[eta0]]],
fpRaw=fDerivs[[All,etaM+eta0+Range[etaP]]],
theImat=IdentityMatroix[numDrvs]},
With[{ckI1=kronI[numDrvs,f0Raw],
ckI2=kronI[numDrvs,fpRaw]},
{f0Raw,fpRaw,ckI1,ckI2}
]]]



spIdentityMatrix[nn_Integer]:=
With[{guts=Map[(#* {1,1})->1&,Range[nn] ]},
SparseArray[guts]]


relRowsSplitUp02[relRows_List,theNumRows_Integer,listLen_Integer]:=
With[{blks=Partition[Range[listLen],theNumRows]},
MapIndexed[Intersection[relRows,#1]-(#2[[1]]-1)*theNumRows&,blks]]


rangeForK[etaK_Integer,theSum_Integer]:=Range[etaK]+theSum
bumpByEtaK[theEtaKs_List,theIdx_Integer]:=(theIdx-1)*theEtaKs;

ofGivenLead[theLeadsList_List,theLead_Integer]:=
Select[theLeadsList,#[[2]]==theLead&]




stochAugB[sigB_?MatrixQ,epsVars_List,tVal_Integer]:=
With[{theNumCols=numCols[sigB],numErrs=Length[epsVars],
epsNow=Transpose[{Through[epsVars[tVal]]}]},
blockMatrix[{
{sigB},
{blockMatrix[{{zeroMatrix[numErrs,theNumCols-1],epsNow}}]},
{blockMatrix[{{zeroMatrix[1,theNumCols-1],IdentityMatrix[1]}}]}
}]]



numCols[theMat_?MatrixQ]:=Length[theMat[[1]]]
numRows[theMat_?MatrixQ]:=Length[theMat]
numRows[xx___]:=Print["badArg numRows",xx];
numCols[xx___]:=Print["badArg numCols",xx];




cmpSSubset02[dim_Integer,theNz_List]:=
With[{numDrvs=numInFacet[dim,2]},
With[{lots=Drop[NestList[MultiIndex`algorithm2,Join[Table[0,{dim-1}],{1}],numDrvs],1]},
With[{theChk=Map[nzAndNotSigma[#,theNz]&,lots]},
Flatten[Position[theChk,True]]]]]
nzAndNotSigma[aList_List,theNz_List]:=
With[{notNz=Complement[Range[Length[aList]],theNz]},
And[aList[[-1]]==0,Max[aList[[notNz]]]==0]]



cmpSSubsetSplit02[aModel_?isAMAModel]:=
With[{bigNeq=getEtaM[aModel]+1},
With[{nonSig=cmpSSubset02[bigNeq,Range[bigNeq-1]]},
{nonSig,Complement[Range[numInFacet[bigNeq,2]],nonSig]}]]


coordinateSplit02[theModel_Symbol,
theLinSys:{{firDrvs_?MatrixQ,secDrvs_?MatrixQ},fMat_?MatrixQ},
theConstant_?VectorQ,assumedMomSubs_List:{}]:=
With[{
sSubset=cmpSSubsetSplit02[theModel],
eta0=getEta0[theModel],
etaP=getEtaP[theModel],
modEta=getEta[theModel]},
With[{},
With[{},
With[{ceSels=Flatten[Map[Range[eta0]+(#-1)*eta0&,
sSubset[[1]]]],
ceRowsSels=Flatten[Map[Range[etaP]+(#-1)*etaP&,
sSubset[[1]]]],
stochSels=Flatten[Map[Range[eta0]+(#-1)*eta0&,
sSubset[[2]]]],
stochRowsSels=Flatten[Map[Range[etaP]+(#-1)*etaP&,
sSubset[[2]]]]},
With[{ceFMat=fMat[[ceRowsSels,ceSels]],
stochFMat=fMat[[stochRowsSels,stochSels]]
},
With[{ceInvStuff=
sFyCmpY0YPInv[
{SparseArray[firDrvs],SparseArray[secDrvs]},modEta,eta0,SparseArray[ceFMat]],
stochInvStuff=
sFyCmpY0YPInv[{firDrvs,secDrvs},modEta,eta0,stochFMat]},
With[{
hooey=
cmpCESolution02[Prepend[ceInvStuff,ceFMat],
theConstant[[ceSels]],theModel]},
With[{gooey=
cmpStochWithCE02[theModel,hooey,ceSels,stochSels,
stochInvStuff,theConstant,theLinSys,
assumedMomSubs]},
{hooey,gooey}]]]]]]]]






eqStochSysCon02[theModel_?isAMAModel,assumedMomSubs_List:{}]:=
With[{sBDerivs=getSBDerivs[theModel,assumedMomSubs],
augmentedSubs=Join[makeSomeMomSubs[getEpsVars[theModel]],assumedMomSubs]},
If[Length[sBDerivs]>1,sBDerivs[[2]],(*already done*)
With[{ignore=eqStochSysCon01[theModel,assumedMomSubs],
igAgain=
Expand[conStochBs[theModel]]},
With[{linSysComponents=
cmpStochLinSys02[theModel,augmentedSubs]},
With[{theConstant=
cmpStochConstant02[theModel,augmentedSubs]},
With[{
splitComps=
coordinateSplit02[theModel,{getFDerivs02[theModel],linSysComponents[[1]]},
theConstant,augmentedSubs]},
theModel /: getSplit02[theModel,assumedMomSubs]=splitComps;
]]]]]]



cmpCESolution03[
{fMat_?MatrixQ,y0InvPart_,interMedPart_?MatrixQ},theConstant_?VectorQ,
theModel_Symbol]:=
With[{
eta0=getEta0[theModel],
etaP=getEtaP[theModel]},
With[{backwardLookingSoln=
(sFyDoKronLUBack[eta0,y0InvPart,theConstant])},
With[{secTerm=
SparseArray[backwardLookingSoln]},
With[{firTerm=
SparseArray[(interMedPart)]},
With[{interMed=
firTerm. secTerm},
With[{bFuncDrvSoln=
backwardLookingSoln-interMed},
With[{matSoln=
Transpose[Partition[
Flatten[bFuncDrvSoln],eta0]],
theFuture=
Transpose[Partition[Flatten[fMat. bFuncDrvSoln],etaP]]},
With[{fullMatSoln=
blockMatrix[{{matSoln},{theFuture}}]},
fullMatSoln
]]]]]]]]






cmpStochWithCE03[theModel_Symbol,ceSoln_?MatrixQ,ceSels_List,stochSels_List,
{stochF0Inv_LinearSolveFunction,stochMidBeg_?MatrixQ},theConstant_?VectorQ,theLinSys_List,
assumedMomSubs_List:{}]:=
With[{etaP=getEtaP[theModel],eta0=getEta0[theModel],etaM=getEtaM[theModel]},
With[{yp=(getFDerivs03[theModel][[1]])[[All,Range[etaP]+etaM+eta0]]},
With[{
ceDiff=kronI[numInFacet[getEtaM[theModel]+1,2],yp][[stochSels]] .
theLinSys[[-1]][[All,ceSels]] . Flatten[Transpose[ceSoln[[Range[eta0]]]]]},
With[{theNoFuture=sFyDoKronLUBack[getEta0[theModel],
stochF0Inv,
Flatten[theConstant[[stochSels]]-ceDiff]]},
Transpose[Partition[Flatten[
(theNoFuture -stochMidBeg . theNoFuture)],eta0]]/.assumedMomSubs]]]]




cmpStochConstant03[theModel_Symbol,assumedMomSubs_List:{}]:=
With[{
augmentedSubs=Join[makeSomeMomSubs[getEpsVars[theModel]],assumedMomSubs],
leadsNeeded=getLeadsNeeded[theModel],
modNz=getNzCols[theModel],
sFDerivs=getFDerivs03[theModel],
theStochBs=getStochConBs[theModel],
etaM=getEtaM[theModel],
eta0=getEta0[theModel],
etaP=getEtaP[theModel],
nlags=getNLags[theModel]
},
With[{
laggedVarsRequired=laggedToKeep[modNz,eta0,nlags],
stateVarDrvs=blockMatrix[{{IdentityMatrix[etaM],zeroMatrix[etaM,1]}}],
bFuncsRequired=stochForMult[leadsNeeded,eta0],
bColsRequired=stochSelectMin[{etaM,eta0,etaP},modNz,nlags]
},
With[{
bAndBCompDrvs=
Flatten[Drop[theStochBs,1][[All,
Length[laggedVarsRequired]+Range[eta0],bColsRequired]],1]},
With[{gFuncDrvs=blockMatrix[{{stateVarDrvs},
{bAndBCompDrvs[[bFuncsRequired,All]]}}]},
(*needs subs or else big time hit for sedmodelopt*)
With[{theGMultVals=Expand[
gMults03[gFuncDrvs,etaM+1,augmentedSubs]]//.augmentedSubs},
With[{theConstant=(-1)*sFDerivs[[2]] . theGMultVals},
Flatten[Transpose[theConstant]]
]]]]]]



stochKronSum03[stochBs_List,{etaM_Integer,neq_Integer,etaP_Integer},
modNz_List,nlags_Integer,lilRel_List,assumedMomSubs_List:{}]:=
With[{supRows02=relRowsSplitUp02[lilRel,neq,Length[stochBs]*neq],
supRows03=relRowsSplitUp03[lilRel,neq,Length[stochBs]*neq]},
With[{ssm= stochSelectMin[{etaM,neq,etaP},modNz,nlags],
ltk= laggedToKeep[modNz,eta0,nlags]},
With[{
bPowsLil=
Map[nullToNullSparseArray[#[[Length[ltk]+Range[neq],Length[ltk]+Range[neq]]]]&,stochBs],
lilStochBs=Map[#[[ssm,ssm]]&,stochBs]},
(* need expand sub or error*)
With[{multListLil=
Map[SparseArray[Expand[gMults02[#,etaM+1,assumedMomSubs]]//.assumedMomSubs]&,
Drop[lilStochBs,1]],
nxtMultListLil=
Map[SparseArray[Expand[gMults03[#,etaM+1,assumedMomSubs]]//.assumedMomSubs]&,
Drop[lilStochBs,1]]},
With[{theINuLil=spIdentityMatrix[Length[multListLil[[1]]]],
theNxtINuLil=spIdentityMatrix[Length[nxtMultListLil[[1]]]]},
With[{},
With[{
specRes=SparseArray[Flatten[
MapIndexed[specificKronSum[#1,#2[[1]],
multListLil,bPowsLil]&,supRows02],1]],
specRes03=SparseArray[Flatten[
MapIndexed[specificKronSum[#1,#2[[1]],
nxtMultListLil,bPowsLil]&,supRows03],1]]
},
{specRes}]]]]]]]



cmpStochLinSys03[theModel_Symbol,augmentedSubs_List]:=
With[{
sFDerivs=getFDerivs03[theModel],
theStochBs=getStochConBs[theModel],
modEta=getEta[theModel],
nlags=getNLags[theModel],
neq=getEta0[theModel],
leadsNeeded=getLeadsNeeded[theModel],
modNz=getNzCols[theModel],
stateVarDim=getEtaM[theModel]+1},
With[{drvsInFacet=numInFacet[stateVarDim,2]},
With[{(*relevantRows=stochKeepers[leadsNeeded,neq,drvsInFacet],*)
lilRel=stochKeepers[leadsNeeded,neq,1]},
With[{kronSums=
 stochKronSum03[Drop[Expand[theStochBs],
-1],modEta,modNz,nlags,lilRel,augmentedSubs]},
With[{fMat=kronSums[[1]][[
reOrderCols[leadsNeeded,drvsInFacet]]]},
Join[{fMat}]]]]]]




relRowsSplitUp03[relRows_List,numRows_Integer,listLen_Integer]:=
With[{blks=Partition[Range[listLen],numRows]},
MapIndexed[Intersection[relRows,#1]-(#2[[1]]-1)*numRows&,blks]]



cmpSSubset03[dim_Integer,theNz_List]:=
With[{numDrvs=numInFacet[dim,2]},
With[{lots=Drop[NestList[MultiIndex`algorithm2,Join[Table[0,{dim-1}],{1}],numDrvs],1]},
With[{theChk=Map[nzAndNotSigma[#,theNz]&,lots]},
Flatten[Position[theChk,True]]]]]
nzAndNotSigma[aList_List,theNz_List]:=
With[{notNz=Complement[Range[Length[aList]],theNz]},
And[aList[[-1]]==0,Max[aList[[notNz]]]==0]]



cmpSSubsetSplit03[aModel_?isAMAModel]:=
With[{bigNeq=getEtaM[aModel]+1},
With[{nonSig=cmpSSubset03[bigNeq,Range[bigNeq-1]]},
{nonSig,Complement[Range[numInFacet[bigNeq,2]],nonSig]}]]


coordinateSplit03[theModel_Symbol,
theLinSys:{{firDrvs_?MatrixQ,secDrvs_?MatrixQ},fMat_?MatrixQ},
theConstant_?VectorQ,assumedMomSubs_List:{}]:=
With[{
sSubset=cmpSSubsetSplit03[theModel],
eta0=getEta0[theModel],
etaP=getEtaP[theModel],
modEta=getEta[theModel]},
With[{},
With[{},
With[{ceSels=Flatten[Map[Range[eta0]+(#-1)*eta0&,
sSubset[[1]]]],
ceRowsSels=Flatten[Map[Range[etaP]+(#-1)*etaP&,
sSubset[[1]]]],
stochSels=Flatten[Map[Range[eta0]+(#-1)*eta0&,
sSubset[[2]]]],
stochRowsSels=Flatten[Map[Range[etaP]+(#-1)*etaP&,
sSubset[[2]]]]},
With[{ceFMat=fMat[[ceRowsSels,ceSels]],
stochFMat=fMat[[stochRowsSels,stochSels]]
},
With[{ceInvStuff=
sFyCmpY0YPInv[
{SparseArray[firDrvs],SparseArray[secDrvs]},modEta,eta0,SparseArray[ceFMat]],
stochInvStuff=
sFyCmpY0YPInv[{firDrvs,secDrvs},modEta,eta0,stochFMat]},
With[{
hooey=
cmpCESolution03[Prepend[ceInvStuff,ceFMat],
theConstant[[ceSels]],theModel]},
With[{gooey=
cmpStochWithCE03[theModel,hooey,ceSels,stochSels,
stochInvStuff,theConstant,theLinSys,
assumedMomSubs]},
{hooey,gooey}]]]]]]]]






eqStochSysCon03[theModel_?isAMAModel,assumedMomSubs_List:{}]:=
With[{sBDerivs=getSBDerivs[theModel,assumedMomSubs],
augmentedSubs=Join[makeSomeMomSubs[getEpsVars[theModel]],assumedMomSubs]},
If[Length[sBDerivs]>1,sBDerivs[[2]],(*already done*)
With[{ignore=eqStochSysCon01[theModel,assumedMomSubs],
igAgain=
Expand[conStochBs[theModel]]},
With[{linSysComponents=
cmpStochLinSys03[theModel,augmentedSubs]},
With[{theConstant=
cmpStochConstant03[theModel,augmentedSubs]},
With[{
splitComps=
coordinateSplit03[theModel,{getFDerivs03[theModel],linSysComponents[[1]]},
theConstant,augmentedSubs]},
theModel /: getSplit03[theModel,assumedMomSubs]=splitComps;
]]]]]]


multExtndMult[aMat_?MatrixQ,pth_?MatrixQ,numExts_Integer]:=
Module[{resMat},
With[{pthRow=Length[pth],extRow=Length[aMat],extCol=Length[pth[[1]]]},
resMat=zeroMatrix[Length[pth]+numExts*extRow,extCol];
resMat[[Range[pthRow],All]]=pth;
With[{aRnge=Range[extRow],pRnge=Range[pthRow]},
Do[resMat[[pthRow+(ii-1)*extRow+aRnge]]=
aMat.resMat[[(ii-1)*extRow+pRnge,All]],{ii,numExts}];
resMat]]]

(*
blendKrons[theMats_List]:=
With[{reGrp=Partition[theMats[[2]],2]},
With[{someKrons=Map[Apply[kron,#]&,reGrp]},
Apply[blockMatrix,{Transpose[{someKrons}]}]]](*probably*)
*)


kron[a_,b_]:=blockMatrix[Outer[Times,a,b]](*probably*)



kronI[nn_Integer,secMat_?MatrixQ]:=SparseArray[kron[
SparseArray[IdentityMatrix[nn]],
(*IdentityMatrix[nn],*)
SparseArray[secMat]]]



mace[exprList_List,limVal_?NumberQ]:=
Max[Abs[Chop[Expand[exprList],limVal]]]

snapShot[expVal_String]:=
With[{theVal=Union[Flatten[Chop[Expand[ToExpression[expVal]]]]]},
Print["snapShot ",expVal,Take[Sort[theVal],Min[25,Length[theVal]]]]];


timeResult[msg_String,funcApplication_]:=
With[{begTime=TimeUsed[],maxMem=MaxMemoryUsed[]/1000000.,
memNow=MemoryInUse[]/1000000.},
With[{res=Evaluate[funcApplication]},
With[{endTime=TimeUsed[]},
Print[msg,{endTime-begTime,memNow,maxMem}];
res]]]
SetAttributes[timeResult,HoldAll]




ps[ss_Integer,nu_List,lambda_List]:=
Module[{},
With[{kRes=kForPs[ss,nu,lambda]},
With[{kVals=getKVals[kRes]},
With[{lVals=Map[getLVals[possibleLForPs[ss,nu,#]]&,kVals]},
Flatten[MapThread[distributeEm,{kVals,lVals}],1]]]]]/;
And[ss>0,ss<=Apply[Plus ,nu]]


distributeEm[xx_,yy_List]:=
Module[{},Map[{xx,#}&,yy]]

noneZero[lList_List]:=FreeQ[Map[theSum,lList],0]


possibleLForPs[ss_Integer,nu_List,aKVal_List]:=
Module[{},
With[{aKSum=Map[(Apply[Plus , #])&, aKVal]},
With[{dd=Length[nu],nn=Apply[Plus , nu]},
With[{possibleLs=If[dd==1,Transpose[{Range[0,nn]}],
Flatten[Apply[Outer,Join[{List},Table[Range[0,nn],{dd}]]],dd-1]]},
With[{possibleTrips=KSubsets[possibleLs,ss]},
With[{smlSet=Select[
If[ss==1,Select[possibleTrips,( (aKSum . #)) ==nu&],
Select[possibleTrips,( (aKSum . #)) ==nu&]],noneZero]},
{dd,nn,possibleLs,possibleTrips,smlSet}]]]]]]/;
And[ss>0,ss<=Apply[Plus ,nu],Length[aKVal]==ss]
getLVals[lRes_List]:=lRes[[5]]



kForPs[ss_Integer,nu_List,lambda_List]:=
Module[{},
With[{dd=Length[nu],mm=Length[lambda],nn=Apply[Plus , nu]},
If[ss===1,
With[{kVals={{lambda}}},
With[{},
{dd,mm,nn,kVals}]],
With[{kPartsRaw=Map[Select[Compositions[#,ss],True&]& ,lambda]},
If[Length[lambda]===1,
With[{kVals=Flatten[kPartsRaw,1]/.xx_Integer->{xx}},
With[{},
{dd,mm,nn,Select[kVals,noneZero]}]],
With[{kPartsOuter=Apply[Outer , Join[{List},kPartsRaw,{1}]]},
With[{kVals=Map[Transpose ,
Flatten[kPartsOuter,mm-1]]},
With[{},
With[{},
With[{},
{dd,mm,nn,Select[kVals,noneZero]}]]]]]]]]]]/;
And[ss>0,ss<=Apply[Plus ,nu],(Apply[Plus , lambda])<=(Apply[Plus , nu])]


getD[kRes_List]:=kRes[[1]]
getM[kRes_List]:=kRes[[2]]
getN[kRes_List]:=kRes[[3]]
getKVals[kRes_List]:=kRes[[4]]


pOrderLst1Smaller[lst1_List,lst2_List]:=
If[(Apply[Plus ,lst1]) <(Apply[Plus ,lst2]),True,
If[(Apply[Plus ,lst1]) >(Apply[Plus ,lst2]),False,
With[{diff=lst1-lst2},
With[{pp1=Flatten[Position[Map[(#<0)&,diff],True]],
pp2=Flatten[Position[Map[(#>0)&,diff],True]]},
If[pp1==={},False,If[First[pp1]<First[pp2],True,False,False]]]]]]


allLambda[mm_Integer,nn_Integer]:=Compositions[nn,mm]

aTerm[nu_List,lambda_List,ss_Integer,{},jj_Integer,gg_List]:=0

aTerm[nu_List,lambda_List,ss_Integer,ps_List,jj_Integer,gg_List]:=
With[{nn=Apply[Plus , nu]},
With[{kVals=ps[[1]],lVals=ps[[2]]},
doExp[doMapDrv[gg,lVals[[jj]]],kVals[[jj]]]/
(doFac[kVals[[jj]]]*((doFac[lVals[[jj]]])^theSum[kVals[[jj]]]))]]/;Length[gg]==Length[ps[[1,1]]]


overJ[nu_List,lambda_List,ss_Integer,ps_List,gg_List]:=
(Apply[Times , Map[aTerm[nu,lambda,ss,ps,#,gg]&,Range[ss]]])

overPs[nu_List,lambda_List,ss_Integer,gg_List]:=
With[{allPs=ps[ss,nu,lambda]},
Apply[Plus , (Map[overJ[nu,lambda,ss,#,gg]&,allPs])]]

overS[nu_List,lambda_List,gg_List]:=
With[{nn=Apply[Plus , nu]},
(Apply[Plus ,(Map[overPs[nu,lambda,#,gg]&,Range[nn]])])]


overLambda[nu_List,ff_,gg_List]:=
With[{nn=theSum[nu],mm=Length[gg]},
doFac[nu]*(Apply[Plus , Flatten[(Map[Function[xx,Map[(doDrv[ff,#]*
overS[nu,#,gg])& ,allLambda[mm,xx]]],
Range[nn]])]])]



vecFdbTermForNu[nu_List,gg_List]:=
With[{nn=theSum[nu],mm=Length[gg]},
doFac[nu]*((Map[Function[xx,Map[(overS[nu,#,gg])& 
,allLambda[mm,xx]]],
Range[nn]]))]


vecFdbGMultForNu[nu_List,gg_List]:=
With[{nn=theSum[nu],mm=Length[gg]},
doFac[nu]*((Map[Function[xx,Map[(overS[nu,#,gg])& 
,allLambda[mm,xx]]],
Range[1,nn]]))]




vecFdbTerm[deg_Integer,gArgs_Integer,gg_List]:=
With[{allNuVals=allLambda[gArgs,deg]},
Map[vecFdbTermForNu[#,gg]&,allNuVals]]

vecFdbGMult[deg_Integer,gArgs_Integer,gg_List]:=
With[{allNuVals=allLambda[gArgs,deg]},
With[{raw=Map[vecFdbGMultForNu[#,gg]&,allNuVals]},
MapThread[toJava,raw]]]

toJava[argPyr___]:=Reverse[Map[Reverse,Transpose[{argPyr}]]]



doMapDrv[gg_List,lVal_List]:=Map[doDrv[#,lVal]&,gg]
doFac[nums_List]:=Apply[Times ,(Map[Factorial,nums])]
theSum[xx_List]:=Apply[Plus , xx]
doExp[bot_List,top_List]:=Module[{},
Apply[Times,
MapThread[myPower,{bot,top}]]]/;Length[bot]==Length[top]
myPower[x_,y_]:=If[x===y===0,1,Power[x,y]]


doDrv[gg_Symbol,drvs_List]:=(Apply[Derivative,drvs])[gg][]
doDrv[gg_,drvs_List]:=With[{offset=theIndex[drvs]},gg[[offset]]]



End[] (* End Private Context *)

EndPackage[]
