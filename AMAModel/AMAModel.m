(* Mathematica Package *)

BeginPackage["AMAModel`",{"ProtectedSymbols`","NumericAMA`","SymbolicAMA`"}]

validateAMAModel::usage = "validateAMAModel  "

generateMathAMAModel::usage = "generateMathAMAModel  "

generateLaTeXAMAModel::usage = "generateLaTeXAMAModel  "

numericLinearizeSystemForOBC::usage="numericLinearizeSystemForOBC[eqns_List]"
symbolicLinearizeSystemForOBC::usage="symbolicLinearizeSystemForOBC[eqns_List]"



equationsToMatrix::usage = "equationsToMatrix  "
(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 
stateVariables::usage = "stateVariables  "

validateAMAModel[fName_String]:=
With[{cmd=StringForm[
"java -cp G:/RES2/mathAMA/AndersonMooreAlgorithm/AndersonMooreAlgorithm/ XsdSchemaSaxValidator AMAModel.xsd ``",
fName]},Run[cmd]]


generateMathAMAModel[fName_String]:=
Module[{arg,ig,eqns},
With[{mathFile=mathOutName[fName]},
With[{cmd=StringForm[
"java   org.apache.xalan.xslt.Process  -IN `` -XSL ``AMAModel2Mma.xsl -OUT ``",fName,Global`$installDir,mathFile]},
Print[cmd,mathFile];Run[cmd];arg=Get[mathFile];
{vars, ig, ig, ig, {ig, eqns}, ig, ig} = Global`AMAModelDefinition[arg];
equationsToMatrix[eqns,vars]]]]

mathOutName[fName_String]:=StringReplace[fName,".xml"->".m"]

generateLaTeXAMAModel[fName_String]:=With[{cmd=
StringForm["java   org.apache.xalan.xslt.Process -IN  ``  -XSL ``AMAModel2LaTeX.xsl -OUT ``",
fName,Global`$installDir,latexOutName[fName]]},Print[cmd];Run[cmd]]


latexOutName[fName_String]:=StringReplace[fName,".xml"->".tex"]


stateVariables[eqns_List]:=With[{vars=Variables[eqns]},
	Union[Cases[vars,xx_[ProtectedSymbols`t+yy_.]->xx,Infinity]]]
lagsLeads[eqns_List]:=With[{vars=Variables[eqns]},
	With[{allEm=Sort[Cases[vars,xx_[ProtectedSymbols`t+yy_.]->yy,Infinity]]},
		{First[allEm],Last[allEm]}]]

fullVec[eqns_List]:=
With[{sVars=stateVariables[eqns],lls=Range @@ lagsLeads[eqns]},
Flatten[Through[sVars[ProtectedSymbols`t+#]]&/@lls]]


fullVec[eqns_List,sVars_List]:=
With[{lls=Range @@ lagsLeads[eqns]},
Flatten[Through[sVars[ProtectedSymbols`t+#]]&/@lls]]

equationsToMatrix[eqns_List]:=With[{fv=fullVec[eqns]},
	Transpose[D[eqns,#]&/@fv]]
	

equationsToMatrix[eqns_List,vars_List]:=With[{fv=fullVec[eqns,vars]},
	Transpose[D[eqns,#]&/@fv]]	



numericLinearizeSystemForOBC[eqns_List]:=
Module[{noCnstr=eqns/.{eps[_][_]->0,eqvdIf[_,xx_,_]->xx},
bmat,phimat,fmat},(*Print[noCnstr];*)
With[{hmat=equationsToMatrix[noCnstr]},(*Print[hmat];*)
{ig,ig,ig,ig,qmat,ig,ig,ig}=numericAMA[hmat,1,1];(*Print[zf,hf];*)
Print["need to generalize to actually compute qmat"];
{hmat,qmat,{bmat,phimat,fmat}=numericComputeBPhiF[hmat,qmat]}
]]
symbolicLinearizeSystemForOBC[eqns_List]:=
Module[{noCnstr=eqns/.{eps[_][_]->0,eqvdIf[_,xx_,_]->xx}},(*Print[noCnstr];*)
With[{hmat=equationsToMatrix[noCnstr]},Print["symbolicLinearize:",hmat];
]]
End[] (* End Private Context *)

EndPackage[]
Print["done reading AMAModel"]
