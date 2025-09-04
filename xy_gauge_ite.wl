(* ::Package:: *)

(* ::Input:: *)
(*X={{0,1},{1,0}};*)
(*Y={{0,-I},{I,0}};*)
(*Z={{1,0},{0,-1}};*)
(*id={{1,0},{0,1}};*)
(**)
(*H:=Subscript[\[CapitalLambda], 1](KroneckerProduct[X,X, id] + KroneckerProduct[Z,Z,id] )+Subscript[\[CapitalLambda], 2](KroneckerProduct[id,X,X] + KroneckerProduct[id,Z,Z] );*)


(* ::Input:: *)
(*Subscript[\[CapitalLambda], 1]=Subscript[\[CapitalLambda], 2]=\[Pi]/4;*)


(* ::Input:: *)
(*exact=Table[F=Tr[MatrixExp[-\[Beta] H]];*)
(*energy=Tr[H . MatrixExp[-\[Beta] H]]/F;*)
(*{\[Beta],N[energy]},{\[Beta],Range[0.,0.15,0.015]*2*5}][[All,2]]*)


(* ::Input:: *)
(*Needs["Wolfram`QuantumFramework`"]*)


(* ::Input:: *)
(*sx:= QuantumOperator["SX", "Label"->"\!\(\*SqrtBox[\(X\)]\)"]*)
(*sz:= QuantumOperator[{"S"}, "Label"->"\!\(\*SqrtBox[\(Z\)]\)"]*)
(*inp :=  PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumState"]["UniformMixture"[3]]*)


(* ::Input:: *)
(*ops = {PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*Subscript[\[Gamma], 1],"XXZ"},{1,2,4}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*Subscript[\[Gamma], 1],"ZZZ"},{1,2,5}], PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*Subscript[\[Gamma], 2],"XXZ"},{2,3,6}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*Subscript[\[Gamma], 2],"ZZZ"},{2,3,7}]};*)


(* ::Input:: *)
(*cxx:=PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumCircuitOperator"][{"+"->{4,5,6,7},ops,ops,ops,ops,ops,sx->{4,5,6,7},"H"->{1,2,3}}, "Parameters"->{Subscript[\[Gamma], 1],Subscript[\[Gamma], 2],\[Beta]}]*)


(* ::Input:: *)
(*czz:=PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumCircuitOperator"][{"+"->{4,5,6,7},ops,ops,ops,ops,ops,sx->{4,5,6,7}},"Parameters"->{Subscript[\[Gamma], 1],Subscript[\[Gamma], 2],\[Beta]}]*)


(* ::Input:: *)
(*Hxx = Table[cxx[<|Subscript[\[Gamma], 1]->2Subscript[\[CapitalLambda], 1],Subscript[\[Gamma], 2]->2Subscript[\[CapitalLambda], 2],\[Beta]->-i|>][inp],{i,0.015,0.15,0.015}];*)
(*Hzz = Table[czz[<|Subscript[\[Gamma], 1]->2Subscript[\[CapitalLambda], 1],Subscript[\[Gamma], 2]->2Subscript[\[CapitalLambda], 2],\[Beta]->-i|>][inp],{i,0.015,0.15,0.015}];*)


(* ::Input:: *)
(*Txx = Table[KeyMap[Function[(1-2*#[[1]])],Hxx[[j]]["Probabilities"]],{j,Length[Hxx]}];*)
(*Tzz = Table[KeyMap[Function[(1-2*#[[1]])],Hzz[[j]]["Probabilities"]],{j,Length[Hzz]}];*)


(* ::Input:: *)
(*approx = Prepend[Table[Total[KeyValueMap[Fxx,Txx[[i]]]]+Total[KeyValueMap[Fzz,Tzz[[i]]]],{i,Length[Txx]}],0.]*)


(* ::Input:: *)
(*shots := 10000*)
(*resx=Values[QuantumMeasurementSimulation[#,PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumMeasurementOperator"]/@{"ZZZZZZZ"},shots]][[1]]&/@Hxx;*)
(*resz=Values[QuantumMeasurementSimulation[#,PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumMeasurementOperator"]/@{"ZZZZZZZ"},shots]][[1]]&/@Hzz;*)


(* ::Input:: *)
(*pos = Range[1,128];*)
(*spins = 1-2*IntegerDigits[pos-1, 2, 7];*)


(* ::Input:: *)
(*Fxx[res_]:=Sum[(Subscript[\[CapitalLambda], 1]*res[[pos]]*spins[[pos]][[1]]*spins[[pos]][[2]]*spins[[pos]][[4]] + Subscript[\[CapitalLambda], 2]*res[[pos]]*spins[[pos]][[2]]*spins[[pos]][[3]]*spins[[pos]][[6]])/shots,{pos,1,128}]*)
(*Fzz[res_]:=Sum[(Subscript[\[CapitalLambda], 1]*res[[pos]]*spins[[pos]][[1]]*spins[[pos]][[2]]*spins[[pos]][[5]] + Subscript[\[CapitalLambda], 2]*res[[pos]]*spins[[pos]][[2]]*spins[[pos]][[3]]*spins[[pos]][[7]])/shots,{pos,1,128}]*)


(* ::Input:: *)
(*approx = Prepend[Map[Fzz,resz]+ Map[Fxx,resx],0.]*)


(* ::Input:: *)
(*exact*)
