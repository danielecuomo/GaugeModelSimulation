(* ::Input:: *)
Needs["Wolfram`QuantumFramework`"]

(* ::Input:: *)
inp :=  QuantumState["UniformMixture"[2]]


(* ::Input:: *)
(*ops={QuantumOperator[{"R",b*e,"XZ"},{1,3}],QuantumOperator[{"R",b*e,"XZ"},{2,4}],QuantumOperator[{"R",\[Beta]*\[Gamma],"ZZZ"},{1,2,5}],QuantumOperator[{"R",b^2*e*g/2,"YZZ"},{1,2,6}],QuantumOperator[{"R",b^2*e*g/2,"ZYZ"},{1,2,7}]};*)


(* ::Input:: *)
(*cz:=PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumCircuitOperator"][{"+"->{3,4,5,6,7},ops,ops,ops,ops,ops, "SX"->{3,4,5,6,7}},"Parameters"->{e,g,b}]*)
(*cx:=PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumCircuitOperator"][{"+"->{3,4,5,6,7},ops,ops,ops,ops,ops,"H"->{1,2}, "SX"->{3,4,5,6,7}},"Parameters"->{e,g,b}]*)


(* ::Input:: *)
(*GAMMA = ETA=Pi/4*)


(* ::Input:: *)
(*pos = Range[1,128];*)
(*spins = 1-2*IntegerDigits[pos-1, 2, 7];*)
(*Fz[res_]:=Sum[GAMMA*res[[pos]]*spins[[pos]][[1]]*spins[[pos]][[2]]*spins[[pos]][[5]]/shots,{pos,1,128}]*)
(*Fx[res_]:=Sum[ETA*res[[pos]]*(spins[[pos]][[3]]*spins[[pos]][[1]] + spins[[pos]][[4]]*spins[[pos]][[2]])/shots,{pos,1,128}]*)


(* ::Input:: *)
(*Hz=Table[cz[<|e->2ETA,g->2GAMMA,b->-j|>][inp],{j,0.01,0.15,0.01}];*)


(* ::Input:: *)
(*Hx=Table[cx[<|e->2ETA,g->2GAMMA,b->-j|>][inp],{j,0.01,0.15,0.01}];*)


(* ::Input:: *)
(*shots := 10000*)
(*resz=Values[QuantumMeasurementSimulation[#,QuantumMeasurementOperator/@{"ZZZZZZZ"},shots]][[1]]&/@Hz;*)
(*resx=Values[QuantumMeasurementSimulation[#,QuantumMeasurementOperator/@{"ZZZZZZZ"},shots]][[1]]&/@Hx;*)


(* ::Input:: *)
(*Prepend[Most[Map[Fz,resz] + Map[Fx,resx]],0.];*)
