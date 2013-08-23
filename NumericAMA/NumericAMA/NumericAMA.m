(* Mathematica Package *)

(* Created by the Wolfram Workbench Jun 17, 2013 *)


BeginPackage["NumericAMA`",{"SymbolicAMA`"}]
(* Exported symbols added here with SymbolName::usage *)  
$zeroTol[]:=10^(-10)


numericAvoidInverse::usage=
"avoidInverse[matA_?MatrixQ,matB_?MatrixQ]"

numericComputeBPhiF::usage=
"numericComputeBPhiF[hmat_?MatrixQ,qmat_?MatrixQ]"
numericComputeB::usage=
"numericComputeB[hmat_?MatrixQ,qmat_?MatrixQ]"

numericComputeZeta::usage=
"numericComputeZeta[qmat_?MatrixQ,upsilon_?MatrixQ,vartheta_?MatrixQ]"

numericObservableStructure::usage=
"numericObservableStructure[ilag_Integer,hmat_?MatrixQ,qmat_?MatrixQ]"


numericRightMostAllZeroQ::usage=
"numericRightMostAllZeroQ[dim_,x_]"

numericShiftRightAndRecord::usage=
"numericShiftRightAndRecord[{auxiliaryConditionsSoFar_,hMatPreShifts_}]"

numericComputeAnnihilator::usage=
"numericComputeAnnihilator[amat_]"

numericAnnihilateRows::usage=
"numericAnnihilateRows[hmat_]"


numericAR::usage=
"numericAR[hmat_]"

numericBiDirectionalAR::usage=
"numericBiDirectionalAR[hmat_]"


numericTransitionMatrix::usage=
"numericTransitionMatrix[transF_]"

numericEliminateInessentialLags::usage=
"numericEliminateInessentialLags[AMatrixVariableListPair_]"


numericMapper::usage=
"numericMapper[bigPi_,bigJ0_]"

numericEvExtend::usage=
"numericEvExtend[evMat_,yMat_,mapper_,transf_]"

numericParticularLam::usage=
"numericParticularLam[es_?MatrixQ,n_Integer,eMap_,transf_]"

numericExtendToBasis::usage=
"numericExtendToBasis[matRows_?MatrixQ]"

numericAsymptoticConstraint::usage=
"numericAsymptoticConstraint[zForward_,zBackward_,gammaForward_]"

numericAMA::usage=
"numericAMA[hmat_?MatrixQ,lags_Integer,infoLag_Integer]"


numericAMAVersion::usage=
"numericAMAVersion[]"
numericUnconditionalCovariance::usage=
"numericUnconditionalCovariance[smat_?MatrixQ,shockCov_?MatrixQ]"

golly::gee=
"golly gee message, has `1` **with** `2` **and** `3` **and** `4` **and** `5`"
Off[golly::gee]
Begin["`Private`"] (* Begin Private Context *) 

numericRightMostAllZeroQ[dim_,x_]:=
        With[{lilvec=Take[x,-dim]},
If[Apply[And,(Map[NumberQ,lilvec])],
(Inner[Times,lilvec,lilvec,Plus]/(dim*dim))<=$zeroTol[],
        If[Apply[And, (Map[Simplify[#] === 0&,lilvec])],True,False,False]]]



numericShiftRightAndRecord[{auxiliaryConditionsSoFar_,hMatPreShifts_}]:=
        If[Apply[Or,Map[Function[x,Apply[And ,Map[#===0&,x]]],hMatPreShifts]],
        Throw[{auxiliaryConditionsSoFar,hMatPreShifts},$zeroRow],
        With[{dim=Length[hMatPreShifts]},
        FoldList[If[numericRightMostAllZeroQ[dim,#2],
                {Append[#1[[1]],Drop[#2,-dim]],
                        Append[#1[[2]],RotateRight[#2,dim]]},
                {#[[1]],Append[#1[[2]],#2]}]&,
                        {auxiliaryConditionsSoFar,{}},hMatPreShifts][[-1]]]]



numericComputeAnnihilator[amat_]:=
With[{adim=Length[amat],
  fa=Flatten[amat]},
  If[Apply[And,(Map[NumberQ, fa])]&&(Max[Abs[fa]] <= $zeroTol[]),
        IdentityMatrix[Length[amat]],
        With[{qr=
 Select[QRDecomposition[amat,Pivoting->True][[1]],
 (Not[Apply[And,(Map[NumberQ,#])]]||(Max[Abs[#]]>0))&]},
   If[Length[qr]==adim,qr(*IdentityMatrix[adim]*),
   Join[qr,Drop[
Apply[Join,numericExtendToBasis[qr]],
Length[qr]]]]]]];




numericAnnihilateRows[hmat_]:=
		With[{dims=Dimensions[hmat]},
        With[{zapper=
                numericComputeAnnihilator[
                subMatrix[hmat,{1,dims[[2]]-dims[[1]]+1},dims[[1]]{1,1}]]},
                If[zapper==={},hmat,zapper . hmat]]]

numericAMA[preHmat_,lags_Integer,ilag_Integer]:=
With[{rescaler=IdentityMatrix[Length[preHmat]]
(*(DiagonalMatrix[(1/Norm /@ preHmat)])*)},
With[{hmat=rescaler. preHmat},
With[{hrows=Length[hmat]},
With[{leads=(Length[hmat[[1]]]/hrows) - (1+lags)},
With[{hcompl=If[lags===0,blockMatrix[{{zeroMatrix[hrows],hmat}}],hmat]},
With[{hcomplr=If[leads===0,blockMatrix[{{hmat,zeroMatrix[hrows]}}],hcompl]},
With[{zfHf=numericAR[hcomplr]},
With[{tm=numericTransitionMatrix[zfHf[[2]]]},
With[{liltmLilJs=If[lags===0&&leads===0,{tm,Range[Length[tm]]},
numericEliminateInessentialLags[{tm,Range[Length[tm]]}]]},
With[{evlsLilevcs=Eigensystem[Transpose[liltmLilJs[[1]]]]},
With[{lrgEvls=Flatten[Position[Abs[#]>1+1.0*^-10& /@evlsLilevcs[[1]],True]]},
expandedevc=expanderFunc[liltmLilJs[[2]],Length[tm]]/@evlsLilevcs[[2,lrgEvls]];
With[{qmat=Join[zfHf[[1]],expandedevc]},
If[And[And[leads==0,Length[qmat] =!= Length[hmat]],Length[qmat] =!= leads*Length[hmat]],
	Throw[{"qmat incorrect number rows",keepStuff,evlsLilevcs[[1]]}]];
With[{bmat=If[lags===0&&leads===0,zeroMatrix[hrows],numericComputeB[hmat,qmat]]},
With[{smat=Inverse[rescaler].If[lags===0&&leads===0,hmat,
numericObservableStructure[ilag,hcomplr,qmat]]},
With[{s0inv=Inverse[If[lags===0&&leads===0,smat,subMatrix[smat,{1,1+If[lags===0,1,lags]*hrows},hrows*{1,1}]]]},
{zfHf,tm,liltmLilJs,evlsLilevcs,qmat,bmat[[Range[Length[hmat]]]],smat,s0inv}]]]]]]]]]]]]]]]



expanderFunc[jList_List,dim_Integer]:=
With[{argMap=Map[If[FreeQ[jList,#1],0,latter[gibbe,Position[jList,#1][[1,1]]]]&,Range[dim]]/.gibbe->Slot[1]},
Function[argMap]/.latter->Part]


numericAR[hmat_]:=
Catch[
FixedPoint[numericShiftRightAndRecord[
        {#[[1]],numericAnnihilateRows[#[[2]]]}]&,
        {{},hmat},Length[hmat[[1]]],
        SameTest->(Length[#1[[1]]] ===Length[#2[[1]]]&)],
        $zeroRow,$rankDeficiency[#1]&]




numericTransitionMatrix[transF_]:=
With[{nr=Length[transF]},
With[{nc=Length[transF[[1]]]-nr},
If[nc == nr, 
-numericAvoidInverse[subMatrix[transF,{1,nc+1},{nr,nr}],
subMatrix[transF,{1,1},{nr,nc}]],
blockMatrix[{{
        blockMatrix[{{zeroMatrix[nc-nr,nr],IdentityMatrix[nc-nr]}}]},
  { blockMatrix[{{-numericAvoidInverse[subMatrix[transF,{1,nc+1},{nr,nr}],
subMatrix[transF,{1,1},{nr,nc}]]}}]}}]]]]


numericEliminateInessentialLags[AMatrixVariableListPair_]:=
    Block[{firstzerocolumn,matrixsize},
    matrixsize=Length[AMatrixVariableListPair[[1]]];
    firstzerocolumn=FindFirstZeroColumn[
            AMatrixVariableListPair[[1]]];
    Return[
        If[IntegerQ[firstzerocolumn],
        numericEliminateInessentialLags[
            List[
            ((AMatrixVariableListPair[[1]])[[
                Drop[Range[1,matrixsize],{firstzerocolumn}],
                Drop[Range[1,matrixsize],{firstzerocolumn}]]]),
            Drop[AMatrixVariableListPair[[2]],{firstzerocolumn}]]],
        AMatrixVariableListPair]]]
FindFirstZeroColumn[AMatrix_]:=Block[{columnsumofabs},
    columnsumofabs=Map[Apply[And,#]&,
        Map[Map[# == 0&,#]& , Transpose[AMatrix]]];
    Return[Apply[Min,Position[columnsumofabs,True]]]]




numericMapper[bigPi_,bigJ0_]:=
Function[evMat,
        (numericAvoidInverse[kron[IdentityMatrix[Length[bigJ0]],evMat] - 
                kron[Transpose[bigJ0],IdentityMatrix[Length[evMat]]],
                kron[Transpose[bigPi],IdentityMatrix[Length[evMat]]]]) ]


numericEvExtend[evMat_,yMat_,mapper_,transf_]:=
With[{xmat=mapper[evMat] . Transpose[{Flatten[Transpose[yMat]]}]},
        blockMatrix[{{Transpose[Partition[Flatten[xmat],
                        Length[yMat](*Length[xmat]/Length[evMat]*)]],yMat}}] .
                    transf]




numericParticularLam[es_?MatrixQ,n_Integer,eMap_,transf_]:=
{Flatten[Append[eMap[{{es[[1,n]]}}] . Transpose[{es[[2,n]]}],es[[2,n]]]]} . transf

(*
numericParticularLam[es_?MatrixQ,n_Integer,eMap_]:=Apply[eMap,Prepend[es[[2,n]],es[[1,n]]]]
*)


numericExtendToBasis[matRows_?MatrixQ]:=
With[{firstQR=QRDecomposition[Transpose[matRows],Pivoting->True]},
With[{bigZ=Select[
firstQR[[1]],
 (Not[Apply[And,(Map[NumberQ,#])]]||(Max[Abs[#]]>0))&]},
With[{secondQR=QRDecomposition[Transpose[bigZ] . bigZ - 
  IdentityMatrix[Length[bigZ[[1]]]],Pivoting->True]},
{bigZ,
Select[
secondQR[[1]],
 (Not[Apply[And,(Map[NumberQ,#])]]||(Max[Abs[#]]>0))&]}]]]



numericComputeBPhiF[hmat_?MatrixQ,qmat_?MatrixQ,leads_Integer]:=
With[{qRows=Length[qmat],qCols=If[qmat==={},0,Length[qmat[[1]]]],hRows=Length[hmat]},
With[{hzero=subMatrix[hmat,{1,qCols-qRows+1},hRows*{1,1}],
hplus=If[leads==0,0,subMatrix[hmat,{1,qCols-qRows+1+hRows},{hRows,qRows}]]},
With[{},
If[False,{{},Inverse[hzero],{"notImplemented"}},
With[{bmats=numericAvoidInverse[
-subMatrix[qmat,{1,qCols-qRows+1},qRows*{1,1}],
subMatrix[qmat,{1,1},{qRows,qCols-qRows}]]},
With[{brmats=subMatrix[bmats,{1,qCols-qRows-hRows+1},{qRows,hRows}]},
With[{phimat=Inverse[hzero + If[leads==0,0,hplus . brmats]]},(*Print["qmat for fmat",{qmat,qRows,hRows}];*)
With[{phihp=-phimat . hplus,
slctrmat=If[leads>1,
blockMatrix[{{zeroMatrix[(leads-1)*hRows,hRows]},{IdentityMatrix[hRows]},
{brmats}}],blockMatrix[{{IdentityMatrix[hRows]},
{brmats}}]]},
With[{busyRes= If[leads==0,0,
Nest[{(*Print["hi",leads,"leads"];*)Append[#[[1]],phihp. subMatrix[#[[2]],{1,1},{(leads)*hRows,hRows}]],
Drop[#[[2]],hRows]}&,{{},slctrmat},leads][[1]]]},
(*Print["howdy",Dimensions[busyRes],busyRes];*)
With[{fmat=If[leads==0,zeroMatrix[hRows,hRows],If[leads>1,
		blockMatrix[{{zeroMatrix[qCols-qRows-hRows,hRows],IdentityMatrix[qCols-qRows-hRows]},busyRes}],busyRes[[1]]]]},
{bmats,phimat,fmat}]]]]]]]]]]

numericComputeB[hmat_?MatrixQ,qmat_?MatrixQ]:=
With[{qRows=Length[qmat],qCols=If[qmat==={},0,Length[qmat[[1]]]],hRows=Length[hmat]},
With[{(*,
hplus=subMatrix[hmat,{1,qCols-qRows+1+hRows},{hRows,qRows}]*)},
With[{lags=(qCols-qRows)/hRows},Off[General::stop];Message[golly::gee,lags,leads,qCols,qRows,hRows];On[General::stop];
If[False,{},
With[{qr=subMatrix[qmat,{1,qCols-qRows+1},qRows*{1,1}],
ql=subMatrix[qmat,{1,1},{qRows,qCols-qRows}]},Off[General::stop];Message[golly::gee,Dimensions[qmat],qr,ql,"",""];On[General::stop];
With[{bmats=numericAvoidInverse[
-qr,ql]},
bmats]]]]]]

numericComputeZeta[qmat_?MatrixQ,upsilon_?MatrixQ,vartheta_?MatrixQ]:=
With[{neq=Length[vartheta]},
With[{nleads=Length[qmat]/neq},
With[{nlags=(If[qmat==={},0,Length[qmat[[1]]]]/neq)-nlags},
With[{qr=subMatrix[qmat,{1,neq*nlags+1},nleads*neq*{1,1}]},
With[{allPow=NestList[# . upsilon&,vartheta . upsilon,nleads-1]},
With[{constituent=blockMatrix[Transpose[{allPow}]]},
qr . constituent]]]]]]


numericObservableStructure[ilag_Integer,hmat_?MatrixQ,qmat_?MatrixQ]:=
With[{hrows=Length[hmat],
qrows=Length[qmat],qcols=If[qmat==={},0,Length[qmat[[1]]]]},
With[{bmat=numericComputeB[hmat,qmat]},
With[{hl=subMatrix[hmat,{1,1},{hrows,qcols-qrows+hrows}]},
With[{hpart=If[ilag<2,hl,blockMatrix[{{zeroMatrix[hrows,hrows*(ilag-1)],hl}}]]},
If[bmat==={},hpart,
With[{hfb=subMatrix[hmat,{1,qcols-qrows+1+hrows},{hrows,qrows}] . bmat},
With[{ltau=Length[bmat[[1]]]},Off[General::stop];Message[golly::gee,ltau,Dimensions[bmat],hrows,qcols,qrows];On[General::stop];
With[{bbar=If[ltau==hrows,bmat,
blockMatrix[{{zeroMatrix[ltau-hrows,hrows],IdentityMatrix[ltau-hrows]},
{subMatrix[bmat,{1,1},{hrows,hrows}],subMatrix[bmat,{1,hrows+1},{hrows,ltau-hrows}]}}]]},
With[{bpart=If[ilag<1,blockMatrix[{{zeroMatrix[hrows],hfb}}],
blockMatrix[{{hfb . MatrixPower[bbar,ilag],zeroMatrix[hrows,hrows*(ilag)]}}]]},
hpart+bpart]]]]]]]]]

numericAvoidInverse[matA_?MatrixQ,matB_?MatrixQ]:=
With[{lu=LUDecomposition[matA]},Transpose[LUBackSubstitution[lu,#]&/@Transpose[matB]]];




numericAMAVersion[]:="$Revision: 2.01 $"

Print["done reading NumericAMA"]
End[] (* End Private Context *)

EndPackage[]
