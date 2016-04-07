(* Mathematica Package *)

BeginPackage["AMAModelDefinition`","ProtectedSymbols`"]

makeParamSubs::usage = "makeParamSubs  "

makeSubbedEqns::usage = "makeSubbedEqns  "

makeShockSubs::usage = "makeShockSubs  "

getEqns::usage = "getEqns  "
(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

makeParamSubs[name_String]:=With[
{params=Global`AMAModelDefinition[name][[3]]},
#[[1]]->#[[2]]&/@params]

makeSubbedEqns[name_String]:=
With[{eqns=Global`AMAModelDefinition[name][[5,2]],
subs=Join[makeParamSubs[name],makeShockSubs[name]]},
(eqns//.subs)/.ssValSubs]

makeShockSubs[name_String]:=With[
{shks=Global`AMAModelDefinition[name][[-2]]},
#[[1]]->0&/@shks]

getEqns[name_String]:=Global`AMAModelDefinition[name][[5,2]]

ssValSubs=(xx_)[t + (_.)] :> ToExpression[ToString[xx]<>"SSVal"]

End[] (* End Private Context *)

EndPackage[]
