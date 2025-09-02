(* ::Package:: *)

(* ::Input:: *)
(*Needs["Wolfram`QuantumFramework`"]*)


(* ::Input:: *)
(*sx:= QuantumOperator["SX", "Label"->"\!\(\*SqrtBox[\(X\)]\)"]*)
(*inp :=  PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumState"]["UniformMixture"[2]]*)


(* ::Input:: *)
(*example:=PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumCircuitOperator"][{"+"->{3,4,5,6,7},PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*\[Eta],"XZ"},{1,3}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*\[Eta],"XZ"},{2,4}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*\[Gamma],"ZZZ"},{1,2,5}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]^2*(\[Eta]*\[Gamma])/2,"YZZ"},{1,2,6}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]^2*(\[Eta]*\[Gamma])/2,"ZYZ"},{1,2,7}], sx->{3,4,5,6,7}},"WireLabels"->{""},"MeasurementWireLabel"->"","Parameters"->{\[Eta],\[Gamma],\[Beta]}]*)
(*example["Diagram"]*)


(* ::Input:: *)
(*ops={PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*\[Eta],"XZ"},{1,3}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*\[Eta],"XZ"},{2,4}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]*\[Gamma],"ZZZ"},{1,2,5}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]^2*(\[Eta]*\[Gamma])/2,"YZZ"},{1,2,6}],PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumOperator"][{"R",\[Beta]^2*(\[Eta]*\[Gamma])/2,"ZYZ"},{1,2,7}]};*)


(* ::Input:: *)
(*cz:=PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumCircuitOperator"][{"+"->{3,4,5,6,7},ops,ops,ops,ops,ops, sx->{3,4,5,6,7}},"Parameters"->{\[Eta],\[Gamma],\[Beta]}]*)
(*cx:=PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumCircuitOperator"][{"+"->{3,4,5,6,7},ops,ops,ops,ops,ops,"H"->{1,2}, sx->{3,4,5,6,7}},"Parameters"->{\[Eta],\[Gamma],\[Beta]}]*)


(* ::Input:: *)
(*\[CapitalGamma]= 0.785398163;*)
(*\[CapitalLambda]= 0.785398163;*)


(* ::Input:: *)
(*pos = Range[1,128];*)
(*spins = 1-2*IntegerDigits[pos-1, 2, 7];*)
(*Fz[res_]:=Sum[\[CapitalGamma]*res[[pos]]*spins[[pos]][[1]]*spins[[pos]][[2]]*spins[[pos]][[5]]/shots,{pos,1,128}]*)
(*Fx[res_]:=Sum[\[CapitalLambda]*res[[pos]]*(spins[[pos]][[3]]*spins[[pos]][[1]] + spins[[pos]][[4]]*spins[[pos]][[2]])/shots,{pos,1,128}]*)


(* ::Input:: *)
(*Hz=Table[cz[<|\[Eta]->2\[CapitalLambda],\[Gamma]->2\[CapitalGamma],\[Beta]->-i|>][inp],{i,0.01,0.15,0.01}];*)


(* ::Input:: *)
(*Hx=Table[cx[<|\[Eta]->2\[CapitalLambda],\[Gamma]->2\[CapitalGamma],\[Beta]->-i|>][inp],{i,0.01,0.15,0.01}];*)


(* ::Input:: *)
(*shots := 10000*)
(*resz=Values[QuantumMeasurementSimulation[#,PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumMeasurementOperator"]/@{"ZZZZZZZ"},shots]][[1]]&/@Hz;*)
(*resx=Values[QuantumMeasurementSimulation[#,PacletSymbol["Wolfram/QuantumFramework", "Wolfram`QuantumFramework`QuantumMeasurementOperator"]/@{"ZZZZZZZ"},shots]][[1]]&/@Hx;*)


(* ::Input:: *)
(*approx = Prepend[Most[Map[Fz,resz] + Map[Fx,resx]],0.]*)


(* ::Input:: *)
(*X={{0,1},{1,0}};*)
(*Z={{1,0},{0,-1}};*)
(*I2={{1,0},{0,1}};*)
(**)
(*H=\[CapitalGamma] KroneckerProduct[Z,Z] +\[CapitalLambda](KroneckerProduct[X,I2]+KroneckerProduct[I2,X]);*)
(**)
(*exact=Most[Table[Z=Tr[MatrixExp[-\[Beta] H]];energy=Tr[H . MatrixExp[-\[Beta] H]]/Z;*)
(*{\[Beta],N[energy]},{\[Beta],Range[0,0.15,0.01]*2*5}][[All,2]]];*)


(* ::Input:: *)
(*gs=Min[Eigenvalues[H]]*)


(* ::Input:: *)
(*p1=Table[{(i-1)*10^-1,approx[[i]]},{i,2,Length[approx]}];*)
(*p2=Table[{(i-1)*10^-1,exact[[i]]},{i,2,Length[exact]}];*)
(*fit=FindFit[p1,yInf+(y0-yInf) Exp[-\[Gamma] x],{{yInf,-1.6},{y0,0},{\[Gamma],0.5}},x];*)
(**)
(*f1[x_]=(yInf+(y0-yInf) Exp[-\[Gamma] x])/. fit;*)
(**)
(*f2=Interpolation[p2,InterpolationOrder->2];*)
(**)
(*xmin=Min[Min[p1[[All,1]]],Min[p2[[All,1]]]];*)
(*xmax=Max[Max[p1[[All,1]]],Max[p2[[All,1]]]];*)
(**)
(*ymin=Min[Min[p1[[All,2]]],Min[p2[[All,2]]]];*)
(*ymax=Max[Max[p1[[All,2]]],Max[p2[[All,2]]]];*)
(**)
(*marginFactor=0.15;*)
(*margin=marginFactor*Abs[gs]+10^-1;*)
(**)
(*ylo=Min[ymin,gs-margin];*)
(*yhi=ymax;*)
(**)
(*Plot[{f1[x],f2[x],gs},{x,xmin,xmax},PlotRange->{{xmin,xmax},{ylo,yhi}},PlotRangePadding->Scaled[.01],PlotStyle->{{Directive[Blue,Dashed,Opacity[0.5]]},{Directive[Red,Dashed,Opacity[0.5]]},{Directive[GrayLevel[0.25],Dashed,Opacity[0.4],AbsoluteThickness[2]]} },PlotLegends->Placed[{"approx","exact","ground"},{0.8,0.65}],AxesLabel->{"\[Beta]","\!\(\*SubscriptBox[\(\:27e8H\:27e9\), \(\[Beta]\)]\)"},GridLines->{Automatic,None},LabelStyle->Directive[FontSize->15],Epilog->{{Blue,Opacity[0.5],PointSize[0.015],Point[p1]},{Red,Opacity[0.5],PointSize[0.015],Point[p2]}}]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
