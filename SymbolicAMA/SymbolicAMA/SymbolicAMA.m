(* Mathematica Package *)

(* Created by the Wolfram Workbench Jun 17, 2013 *)

(* Mathematica Package *)

BeginPackage["SymbolicAMA`"]
(* Exported symbols added here with SymbolName::usage *)  
$rankDeficiency::usage="symbol for error return"
vec::usage=
"vec[mat_List]"





computeXNext::usage="computeXNext[xInit_?MatrixQ,bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]"

computeXPath::usage="computeXPath[xInit_?MatrixQ,bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]"

computeDelXPath::usage="computeXPath[bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]"

zeroMatrix::usage="zeroMatrix[dims]"
subMatrix::usage="subMatrix[hmat,offset,rowscols]"
blockMatrix::usage="blockMatrix[{{aa,bb,cc}}]"
symbolicComputeVartheta::usage=
"symbolicComputeVartheta[phimat_List,psimat_List,fmat_List,upsilon_List]"

symbolicUnconditionalCovariance::usage=
"symbolicUnconditionalCovariance[smat_List,shockCov_List]"

avoidInverse::usage=
"avoidInverse[matA_List,matB_List]"

symbolicComputeBPhiF::usage=
"symbolicComputeBPhiF[hmat_List,qmat_List]"
symbolicComputeB::usage=
"symbolicComputeB[hmat_List,qmat_List]"


symbolicObservableStructure::usage=
"symbolicObservableStructure[ilag_Integer,hmat_List,qmat_List]"

symbolicRightMostAllZeroQ::usage=
"symbolicRightMostAllZeroQ[dim_,x_]"


symbolicShiftRightAndRecord::usage="symbolicShiftRightAndRecord[{auxiliaryConditionsSoFar_,hMatPreShifts_}]
"
symbolicAnnihilateRows::usage="symbolicAnnihilateRows[hmat_]"

symbolicComputeAnnihilator::usage="symbolicComputeAnnihilator[amat_]"

symbolicBiDirectionalAR::usage="symbolicBiDirectionalAR[hmat_]"

symbolicAR::usage="symbolicAR[hmat_]"

symbolicEliminateInessentialLags::usage="symbolicEliminateInessentialLags[AMatrixVariableListPair_]"

symbolicSqueeze::usage="symbolicSqueeze[auxFAuxB_,transMat_]"

symbolicAvoidInverse::usage="symbolicAvoidInverse[matA_List,matB_List]"

symbolicMapper::usage="symbolicMapper[bigPi_,bigJ0_]"

symbolicEvExtend::usage="symbolicEvExtend[evMat_,yMat_,mapper_,transf_]"

symbolicParticularLam::usage="symbolicParticularLam[es_List,n_Integer,eMap_,transf_]"

symbolicExtendToBasis::usage="symbolicExtendToBasis[matRows_List]"

symbolicTransitionMatrix::usage="symbolicTransitionMatrix[transF_]"

kron::usage="kron[a_,b_]"

symbolicAMAVersion::usage="symbolicAMAVersion[]"

symbolicEigenvecSeries::usage = "symbEigenvecSeries  "

symbolicEigNS::usage = "symbolicEigNS  "

seriesAVal::usage = "seriesAVal  "

genericQmat::usage="genericQmat"

sparseToDenseMat::usage="sparseToDenseMat"

denseColToSparseMat::usage="denseColToSparseMat"

denseToSparseMat::usage="denseToSparseMat"


Begin["`Private`"] (* Begin Private Context *) 

zeroMatrix[anInt_Integer]:=ConstantArray[0,{anInt,anInt}]
zeroMatrix[intA_Integer,intB_Integer]:=ConstantArray[0,{intA,intB}]

blockMatrix[aMat_]:=ArrayFlatten[aMat]

subMatrix[mat_?MatrixQ,offset_List,dims_List]:=
With[{theRows={offset[[1]],offset[[1]]+dims[[1]]-1},
theCols={offset[[2]],offset[[2]]+dims[[2]]-1}},
Take[mat,theRows,theCols]]

symbolicRightMostAllZeroQ[dim_,x_]:=
With[{lilvec=Take[x,-dim]},
        If[Apply[And, (Map[Simplify[#] == 0&, lilvec])],
                True,False,False]]

genericQmat[needed_Integer,
paramSubs_List,evls_?VectorQ,evcs_?MatrixQ,zf_?MatrixQ]:=
If[evcs===$Failed,$Failed,
With[{bigEvals=Flatten[Position[Abs[#]>1& /@(evls/.paramSubs),True]]},
If[Length[bigEvals]=!=needed-Length[zf],$Failed,Join[zf,evcs[[bigEvals]]]]]]



genericQmat[___,$Failed,___]:=$Failed


(*
symbolicRightMostAllZeroQ[dim_,x_]:=
With[{lilvec=Take[x,-dim]},
        If[Chop[Norm[FullSimplify[lilvec],2]] == 0,
                True,False,False]]
*)

symbolicShiftRightAndRecord[{auxiliaryConditionsSoFar_,hMatPreShifts_}]:=
        If[Apply[Or,Map[Function[x,Apply[And ,Map[#===0&,x]]],hMatPreShifts]],
        Throw[{auxiliaryConditionsSoFar,hMatPreShifts},$zeroRow],
        With[{dim=Length[hMatPreShifts]},
        FoldList[If[symbolicRightMostAllZeroQ[dim,#2],
                {Append[#1[[1]],Drop[#2,-dim]],
                        Append[#1[[2]],RotateRight[#2,dim]]},
                {#[[1]],Append[#1[[2]],#2]}]&,
                        {auxiliaryConditionsSoFar,{}},hMatPreShifts][[-1]]]]




symbolicAnnihilateRows[hmat_]:=With[{dims=Dimensions[hmat]},
        With[{zapper=
                symbolicComputeAnnihilator[
                subMatrix[hmat,{1,dims[[2]]-dims[[1]]+1},dims[[1]]{1,1}]]},
                If[zapper==={},hmat,zapper . hmat]]]



symbolicComputeAnnihilator[amat_]:=
With[{dim=Length[amat]},If[MatrixRank[amat]<Length[amat],
        With[{ns=RowReduce[blockMatrix[{{amat,IdentityMatrix[dim]}}]]},
        subMatrix[ns,{1,dim+1},{dim,dim}]],{}]]


symbolicBiDirectionalAR[hmat_]:=
{symbolicAR[hmat],Map[(Map[Reverse,#])& ,symbolicAR[Map[Reverse,hmat]]]};


symbolicAR[hmat_]:=
Catch[
FixedPoint[symbolicShiftRightAndRecord[
        {#[[1]],Simplify[symbolicAnnihilateRows[#[[2]]]]}]&,
        {{},hmat},Length[hmat[[1]]],
        SameTest->(Length[#1[[1]]] ===Length[#2[[1]]]&)],
        $zeroRow,$rankDeficiency[#1]&]
(*

expandedEVecs=Array[0&,{1,4}]
expandedEVecs[[All,positions]]=lileVecs[[{2}]]


*)


Unprotect[TimeConstrained];(*TimeConstrained[$Failed,___]:=$Failed;*)

TimeConstrained[
  symbolicComputeBPhiF[_, $Failed], ___] := {$Failed, $Failed, \
$Failed}; Protect[TimeConstrained];

symbolicEliminateInessentialLags[AMatrixVariableListPair_]:=
    Block[{firstzerocolumn,matrixsize},
    matrixsize=Length[AMatrixVariableListPair[[1]]];
    firstzerocolumn=FindFirstZeroColumn[
            AMatrixVariableListPair[[1]]];
    Return[
        If[IntegerQ[firstzerocolumn],
        symbolicEliminateInessentialLags[
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






symbolicSqueeze[auxFAuxB_,transMat_]:=
Module[{},
With[{extend=symbolicExtendToBasis[auxFAuxB]},(*{{p1,p2},mat1,mat2}*)
With[{r12= extend[[2]] . Transpose[extend[[1,2]]],
zeroDim=Length[extend[[1,1]]],
nonZeroDim=Length[extend[[1,2]]],
pmat=Join[extend[[1,1]],extend[[1,2]]]},
With[{ptrans=pmat .transMat. Transpose[pmat]},
{pmat ,subMatrix[ptrans,{1,1},zeroDim*{1,1}],
subMatrix[ptrans,{(zeroDim+1),1},{nonZeroDim,zeroDim}],
subMatrix[ptrans,(zeroDim+1){1,1},nonZeroDim*{1,1}] -
subMatrix[ptrans,{(zeroDim+1),1},{nonZeroDim,zeroDim}]  . r12 }]]]]

symbolicAvoidInverse[matA_List,matB_List]:=Inverse[matA] . matB
(*
With[{lu=LUDecomposition[matA]},Transpose[LUBackSubstitution[lu,#]&/@Transpose[matB]]];
*)

symbolicUnconditionalCovariance[smat_List,shockCov_List]:=
With[{srows=Length[smat],scols=Length[smat[[1]]]},
With[{sinv=LUDecomposition[subMatrix[smat,{1,scols-srows+1},{srows,srows}]]},
If[scols - 2 * srows >0,
With[{psiMat=blockMatrix[{{zeroMatrix[scols-(2* srows),scols-srows]},
{zeroMatrix[srows,scols-2*srows],LUBackSubstitution[sinv,
 Transpose[LUBackSubstitution[sinv , shockCov]]]}}]},Print["got psiMat(true)"];
With[{varTrans=
  blockMatrix[{{zeroMatrix[scols-2*srows,srows],IdentityMatrix[scols-2*srows]},
    {-LUBackSubstitution[sinv , subMatrix[smat,{1,1},{srows,scols-srows}]]}}]},Print["got varTrans(true)"];
Partition[
LUBackSubstitution[
LUDecomposition[(IdentityMatrix[(scols-srows)^2]-kron[varTrans,varTrans])],Flatten[psiMat]],scols-srows]]],
With[{psiMat=LUBackSubstitution[sinv,
  Transpose[LUBackSubstitution[sinv ,shockCov]]]},Print["got psiMat(false)"];
With[{varTrans= -LUBackSubstitution[sinv , subMatrix[smat,{1,1},{srows,scols-srows}]]},Print["got varTrans(false)"];
With[{guts=(IdentityMatrix[(scols-srows)^2]-kron[varTrans,varTrans])},Print["got guts"];
With[{decomp=LUDecomposition[guts]},Print["got decomp"];
{Partition[
LUBackSubstitution[decomp,
Flatten[psiMat]],scols-srows],guts,Flatten[psiMat]}]]]]]]]



symbolicMapper[bigPi_,bigJ0_]:=
Function[evMat,
        symbolicAvoidInverse[kron[IdentityMatrix[Length[bigJ0]],evMat] - 
                kron[Transpose[bigJ0],IdentityMatrix[Length[evMat]]] ,
                kron[Transpose[bigPi],IdentityMatrix[Length[evMat]]]]]


symbolicEvExtend[evMat_,yMat_,mapper_,transf_]:=
With[{xmat=mapper[evMat] . Transpose[{Flatten[Transpose[yMat]]}]},
        blockMatrix[{{Transpose[Partition[Flatten[xmat],
                        Length[yMat](*Length[xmat]/Length[evMat]*)]],yMat}}].
                    transf]




symbolicParticularLam[es_List,n_Integer,eMap_,transf_]:=
{Flatten[Append[eMap[{{es[[1,n]]}}] . Transpose[{es[[2,n]]}],es[[2,n]]]]} . transf



symbolicExtendToBasis[matRows_?MatrixQ]:=
With[{spaceDim=Length[matRows[[1]]]},
With[{idim=IdentityMatrix[spaceDim]},
With[{rowEsch=Select[RowReduce[matRows],!symbolicRightMostAllZeroQ[spaceDim,#]&]},
With[{cCols=Map[Min[Flatten[Position[#,1][[1]]]]& , rowEsch]},
With[{notCCols=Complement[Range[spaceDim],cCols]},
        {{idim[[cCols]],
idim[[notCCols]]},rowEsch,idim[[notCCols]]}]]]]]/;matRows !={}


symbolicTransitionMatrix[transF_]:=
With[{nr=Length[transF]},
With[{nc=Length[transF[[1]]]-nr},
If[nc == nr, 
-symbolicAvoidInverse[subMatrix[transF,{1,nc+1},{nr,nr}] ,
subMatrix[transF,{1,1},{nr,nc}]],
With[{topPart=blockMatrix[{{zeroMatrix[nc-nr,nr],IdentityMatrix[nc-nr]}}],
botPart=blockMatrix[{{-symbolicAvoidInverse[subMatrix[transF,{1,nc+1},{nr,nr}],
subMatrix[transF,{1,1},{nr,nc}]]}}]},
blockMatrix[{{topPart},{botPart}}]]]]]


kron[a_,b_]:=blockMatrix[Outer[Times,a,b]]








symbolicComputeVartheta[phimat_?MatrixQ,psimat_?MatrixQ,fmat_?MatrixQ,upsilon_?MatrixQ]:=
With[{fRows=Length[fmat],hRows=Length[phimat],zLen=Length[psimat[[1]]],
     krnprt=kron[Transpose[upsilon],fmat]},
With[{guts=
If[fRows>hRows,
avoidInverse[IdentityMatrix[Length[krnprt]] - krnprt,
vec[blockMatrix[{{zeroMatrix[fRows-hRows,zLen]},{phimat. psimat}}]]],
avoidInverse[IdentityMatrix[Length[krnprt]] - krnprt,vec[phimat. psimat]]]},
(Transpose[Partition[Flatten[guts],Length[fmat]]])[[Range[fRows-hRows+1,fRows]]]
]]

avoidInverse[matA_?MatrixQ,matB_?MatrixQ]:=
With[{lu=LUDecomposition[matA]},Transpose[LUBackSubstitution[lu,#]&/@Transpose[matB]]];


vec[mat_?MatrixQ]:=Transpose[{Flatten[Transpose[mat]]}];

symbolicEigenvecSeries[aMat_?MatrixQ,lam_,varDegsPos_?MatrixQ]:=
With[{},
	With[{nSpace=symbolicEigNS[aMat,lam,varDegsPos]},
Map[seriesAVal[#,varDegsPos]&,nSpace,{2}]]]


symbolicEigNS[aMat_?MatrixQ,lam_,varDegsPos_?MatrixQ]:=
With[{dim=Dimensions[aMat][[1]]},
	With[{nSpace=NullSpace[Transpose[aMat]-lam*IdentityMatrix[dim]]},nSpace]]

seriesAVal[theVal_,varDegsPos_?MatrixQ]:=
Series@@ Prepend[varDegsPos,theVal]

symbolicComputeB[___,$Failed,___]:=$Failed
symbolicComputeBPhiF[___,$Failed,___]:={$Failed,$Failed,$Failed}



symbolicComputeBPhiF[hmat_?MatrixQ,qmat_?MatrixQ]:=
With[{qRows=Length[qmat],qCols=Length[qmat[[1]]],hRows=Length[hmat]},
With[{leads=qRows/hRows,hzero=subMatrix[hmat,{1,qCols-qRows+1},hRows*{1,1}],
hplus=subMatrix[hmat,{1,qCols-qRows+1+hRows},{hRows,qRows}]},
With[{bmats=avoidInverse[
-subMatrix[qmat,{1,qCols-qRows+1},qRows*{1,1}],
subMatrix[qmat,{1,1},{qRows,qCols-qRows}]]},
With[{brmats=subMatrix[bmats,{1,qCols-qRows-hRows+1},{qRows,hRows}]},
With[{phimat=Inverse[hzero + hplus . brmats]},
With[{phihp=-phimat . hplus,
slctrmat=If[leads>1,
blockMatrix[{{zeroMatrix[(leads-1)*hRows,hRows]},{IdentityMatrix[hRows]},
{brmats}}],blockMatrix[{{IdentityMatrix[hRows]},
{brmats}}]]},
With[{busyRes= 
Nest[{(*Print["hi",leads,"leads"];*)Append[#[[1]],phihp. subMatrix[#[[2]],{1,1},{(leads)*hRows,hRows}]],
Drop[#[[2]],hRows]}&,{{},slctrmat},leads][[1]]},(*Print["howdy",Dimensions[busyRes],busyRes];*)
With[{fmat=If[leads>1,
With[{theBottom=blockMatrix[{busyRes}]},
blockMatrix[{{
		blockMatrix[{{zeroMatrix[qRows-hRows,hRows],IdentityMatrix[qRows-hRows]}}]},{theBottom}}]],busyRes[[1]]]},
{bmats,phimat,fmat}]]]]]]]]

symbolicComputeB[hmat_?MatrixQ,qmat_?MatrixQ]:=
With[{qRows=Length[qmat],qCols=Length[qmat[[1]]]},
With[{},
With[{bmats=avoidInverse[
-subMatrix[qmat,{1,qCols-qRows+1},qRows*{1,1}],
subMatrix[qmat,{1,1},{qRows,qCols-qRows}]]},
bmats]]]



symbolicObservableStructure[ilag_Integer,hmat_?MatrixQ,qmat_?MatrixQ]:=
With[{hrows=Length[hmat],
qrows=Length[qmat],qcols=Length[qmat[[1]]]},
With[{bmat=symbolicComputeB[hmat,qmat]},
With[{hfb=subMatrix[hmat,{1,qcols-qrows+1+hrows},{hrows,qrows}] . bmat},
With[{ltau=Length[bmat[[1]]]},
With[{bbar=If[ltau==hrows,bmat,
blockMatrix[{{zeroMatrix[ltau-hrows,hrows],IdentityMatrix[ltau-hrows]},
{subMatrix[bmat,{1,1},{hrows,ltau}]}}]]},
With[{hl=subMatrix[hmat,{1,1},{hrows,qcols-qrows+hrows}]},
With[{hpart=If[ilag<2,hl,blockMatrix[{{zeroMatrix[hrows,hrows*(ilag-1)],hl}}]],
bpart=If[ilag<1,blockMatrix[{{zeroMatrix[hrows],hfb}}],
blockMatrix[{{hfb . MatrixPower[bbar,ilag],zeroMatrix[hrows,hrows*(ilag)]}}]]},
hpart+bpart]]]]]]]


symbolicObservableStructureBOneLagLead[ilag_Integer,hmat_?MatrixQ,bmat_?MatrixQ]:=
With[{hrows=Length[hmat]},
With[{},
With[{hfb=subMatrix[hmat,{1,hrows+1+hrows},{hrows,hrows}] . bmat},
With[{ltau=Length[bmat[[1]]]},
With[{bbar=If[ltau==hrows,bmat,
blockMatrix[{{zeroMatrix[ltau-hrows,hrows],IdentityMatrix[ltau-hrows]},
{subMatrix[bmat,{1,1},{hrows,ltau}]}}]]},
With[{hl=subMatrix[hmat,{1,1},{hrows,hrows+hrows}]},
With[{hpart=If[ilag<2,hl,blockMatrix[{{zeroMatrix[hrows,hrows*(ilag-1)],hl}}]],
bpart=If[ilag<1,blockMatrix[{{zeroMatrix[hrows],hfb}}],
blockMatrix[{{hfb . MatrixPower[bbar,ilag],zeroMatrix[hrows,hrows*(ilag)]}}]]},
hpart+bpart]]]]]]]

symbolicAMAVersion[]:="$Revision: 2.0 $ $Date: 2011/02/04 $"




computeDelXPath[bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]:=
With[{neq=Length[bMat],numz=Length[psiMat[[1]]]},
With[{nleads=Length[fMat]/neq,numNonZeroZ=Length[zPath]/numz},
With[{preMat=makePreMat[neq,nleads],
postMats=makePostMat[phiMat,psiMat,#,nleads]&/@Partition[zPath,numz],
fPows=NestList[fMat .#&,IdentityMatrix[Length[fMat]],numNonZeroZ-1]},
With[{theProds=MapThread[#1 .#2&,{fPows,postMats}]},
preMat . (Plus @@ theProds)]]]]/;
And[isSquare[phiMat],isSquare[fMat],Length[bMat]==Length[phiMat],
Mod[Length[fMat],Length[phiMat]]==0,Length[psiMat]==Length[phiMat],
Mod[Length[zPath],Length[psiMat[[1]]]]==0]


computeXNext[xInit_?MatrixQ,bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]:=
With[{},
With[{
theDel=computeDelXPath[bMat,phiMat,fMat,psiMat,zPath]},
bMat . xInit + theDel]]


computeXPath[xInit_?MatrixQ,pathLength_Integer,bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]:=
With[{},
With[{
theDel=computeDelXPath[bMat,phiMat,fMat,psiMat,zPath]},
bMat . xInit + theDel]]

computeXPath[xInit_?MatrixQ,0,bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]:={}

computeXPath[xInit_?MatrixQ,pathLength_Integer,bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]:=
With[{nxtX=computeXNext[xInit,bMat,phiMat,fMat,psiMat,zPath]},
Join[nxtX,computeXPath[updateXinit[xInit,nxtX],
pathLength-1,bMat,phiMat,fMat,psiMat,updateZ[zPath,Length[psiMat[[1]]]]]]]

updateXinit[oldX_?MatrixQ,nxtX_?MatrixQ]:=
Join[Drop[oldX,Length[nxtX]],nxtX]


updateZ[oldZ_?MatrixQ,numZ_Integer]:=
Join[Drop[oldZ,numZ],ConstantArray[0,{numZ,1}]]





chkComputeDelXPath[bMat_?MatrixQ,phiMat_?MatrixQ,fMat_?MatrixQ,psiMat_?MatrixQ,zPath_?MatrixQ]:=
And[isSquare[phiMat],isSquare[fMat],Length[bMat]==Length[phiMat],
Mod[Length[fMat],Length[phiMat]]==0,Length[psiMat]==Length[phiMat],
Mod[Length[zPath],Length[psiMat[[1]]]]==0]


vectorOfMats[vom_List]:=
With[{matsQ=MatrixQ/@vom},
If[And @@ matsQ,Length[Union[Dimensions /@vom]]==1,False]]




makePreMat[neq_Integer,1]:=IdentityMatrix[neq]
makePreMat[neq_Integer,nleads_Integer]:=
blockMatrix[{{ConstantArray[0,{neq,neq*(nleads-1)}],IdentityMatrix[neq]}}]

makePostMat[phiMat_?MatrixQ,psiMat_?MatrixQ,zMat_?MatrixQ,nleads_Integer]:=
With[{neq=Length[phiMat],zCols=Length[zMat[[1]]]},
blockMatrix[{{ConstantArray[0,{neq*(nleads-1),zCols}]},{phiMat.psiMat.zMat}}]]

makePostMat[phiMat_?MatrixQ,psiMat_?MatrixQ,zMat_?MatrixQ,1]:=
With[{},
phiMat.psiMat.zMat]




isSquare[mat_?MatrixQ]:=Length[mat]==Length[mat[[1]]]
isSquare[{{}}]=False


presprow[aRow_List]:=
With[{allnz=Complement[Range[Length[aRow]],Flatten[Position[Chop[aRow],0]]]},
With[{allv=aRow[[allnz]]},{allv,allnz}]]


denseToSparseMat[mat_List]:=
With[{res=presprow/@ mat},
With[{lens=Length[#[[1]]]&/@res},Append[blockMatrix @@ {{res}},1+FoldList[Plus,0,lens]]]]


denseColToSparseMat[vec_List]:=
{vec,Table[1,{Length[vec]}],Range[Length[vec]+1]}

sparseToDenseMat[rows_,cols_,a_List,ja_List,ia_List]:=
Module[{s=zeroMatrix[rows,cols]},
With[{rowBounds=Range[#[[1]],#[[2]]-1]&/@Partition[ia,2,1]},
With[{assn=Transpose /@ 
		Transpose[{Part[ja,#1] & /@ rowBounds,Part[a,#1] & /@ rowBounds}]},
MapIndexed[Function[{x,pos},
Map[(s[[pos[[1]],#1[[1]]]]=#1[[2]])&,x]],assn];
s
]]]

End[] (* End Private Context *)

EndPackage[]
Print["done reading SymbolicAMA"]
