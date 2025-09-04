Needs["Wolfram`QuantumFramework`"]

inp :=  QuantumState["UniformMixture"[3]]

ops = {QuantumOperator[{"R",\[Beta]*g1,"XXZ"},{1,2,4}],QuantumOperator[{"R",b*g1,"ZZZ"},{1,2,5}], QuantumOperator[{"R",b*g1,"XXZ"},{2,3,6}],QuantumOperator[{"R",b*g2,"ZZZ"},{2,3,7}]};

cxx:=QuantumCircuitOperator[{"+"->{4,5,6,7},ops,ops,ops,ops,ops,"SX"->{4,5,6,7},"H"->{1,2,3}}, "Parameters"->{g1,g2,b}]

czz:=QuantumCircuitOperator[{"+"->{4,5,6,7},ops,ops,ops,ops,ops,"SX"->{4,5,6,7}},"Parameters"->{g1,g2,b}]

Hxx = Table[cxx[<|g1->2G1,g2->2G2,b->-j|>][inp],{j,0.015,0.15,0.015}];
Hzz = Table[czz[<|g1->2G1,g2->2G2,b->-j|>][inp],{j,0.015,0.15,0.015}];

Txx = Table[KeyMap[Function[(1-2*#[[1]])],Hxx[[j]]["Probabilities"]],{j,Length[Hxx]}];
Tzz = Table[KeyMap[Function[(1-2*#[[1]])],Hzz[[j]]["Probabilities"]],{j,Length[Hzz]}];

Prepend[Table[Total[KeyValueMap[Fxx,Txx[[j]]]]+Total[KeyValueMap[Fzz,Tzz[[j]]]],{j,Length[Txx]}],0.];

shots := 10000
resx=Values[QuantumMeasurementSimulation[#,QuantumMeasurementOperator/@{"ZZZZZZZ"},shots]][[1]]&/@Hxx;
resz=Values[QuantumMeasurementSimulation[#,QuantumMeasurementOperator/@{"ZZZZZZZ"},shots]][[1]]&/@Hzz;

pos = Range[1,128];
spins = 1-2*IntegerDigits[pos-1, 2, 7];

Fxx[res_]:=Sum[(G1*res[[pos]]*spins[[pos]][[1]]*spins[[pos]][[2]]*spins[[pos]][[4]] + G2*res[[pos]]*spins[[pos]][[2]]*spins[[pos]][[3]]*spins[[pos]][[6]])/shots,{pos,1,128}]
Fzz[res_]:=Sum[(G1*res[[pos]]*spins[[pos]][[1]]*spins[[pos]][[2]]*spins[[pos]][[5]] + G2*res[[pos]]*spins[[pos]][[2]]*spins[[pos]][[3]]*spins[[pos]][[7]])/shots,{pos,1,128}]

Prepend[Map[Fzz,resz]+ Map[Fxx,resx],0.];
