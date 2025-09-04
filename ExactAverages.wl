(* ::Input:: *)
X={{0,1},{1,0}};
Y={{0,I},{0,-I}};
Z={{1,0},{0,-1}};
id={{1,0},{0,1}};
(* ::Input:: *)
GAMMA1=GAMMA2=Pi/4;
Hxy:=GAMMA1(KroneckerProduct[X,X, id] + KroneckerProduct[Y,Y,id] )+GAMMA2(KroneckerProduct[id,X,X] + KroneckerProduct[id,Y,Y] );
(* ::Input:: *)
AVGxy=Table[F=Tr[MatrixExp[-b Hxy]];energy=Tr[Hxy . MatrixExp[-b Hxy]]/F;{b,N[energy]},{b,Range[0.,0.15,0.015]*2*5}][[All,2]];
(* ::Input:: *)
GAMMA = ETA = Pi/4;
Hxz=GAMMA*KroneckerProduct[Z,Z] +ETA(KroneckerProduct[X,id]+KroneckerProduct[id,X]);
(* ::Input:: *)
AVGxz=Table[F=Tr[MatrixExp[-b Hxz]];energy=Tr[Hxz . MatrixExp[-b Hxz]]/F;{b,N[energy]},{b,Range[0.,0.15,0.015]*2*5}][[All,2]];
