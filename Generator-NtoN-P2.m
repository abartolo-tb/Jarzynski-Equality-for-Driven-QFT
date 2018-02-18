(* ::Package:: *)

Print["Start Time: "  DateString[Now]]
temperature=1;
Omega[q_]:=Sqrt[1+q^2]
myDiscriminant[W_,q1_,q2_,t1_,t2_,phi_]:= ((q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2])-(W+Omega[q1]+Omega[q2])^2)/(2Max[Abs[(W+Omega[q1]+Omega[q2])],10^-100]))^2+((q1 Cos[t1] + q2 Cos[t2])/Max[Abs[(W+Omega[q1]+Omega[q2])],10^-100])^2-1
DeltaRoots[W_,\[Beta]_,q1_,q2_,t1_,t2_,phi_]:=If[myDiscriminant[W,q1,q2,t1,t2,phi]<= 0,{10^-100},{1/(1-((q1 Cos[t1] + q2 Cos[t2])/Max[Abs[(W+Omega[q1]+Omega[q2])],10^-100])^2) (1/2 (1-(q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2]))/(Max[Abs[(W+Omega[q1]+Omega[q2])],10^-100])^2)(q1 Cos[t1] + q2 Cos[t2])+Sqrt[myDiscriminant[W,q1,q2,t1,t2,phi]]),1/(1-((q1 Cos[t1] + q2 Cos[t2])/Max[Abs[(W+Omega[q1]+Omega[q2])],10^-100])^2) (1/2 (1-(q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2]))/(Max[Abs[(W+Omega[q1]+Omega[q2])],10^-100])^2)(q1 Cos[t1] + q2 Cos[t2])-Sqrt[myDiscriminant[W,q1,q2,t1,t2,phi]])},{10^-100}]
myIntegrand[W_,\[Beta]_,q1_,q2_,q3_,t1_,t2_,phi_]:=Boole[q3>10^-50]*Boole[ Abs[W+Sqrt[1+q1^2]+Sqrt[1+q2^2]-Sqrt[1+q3^2]-Sqrt[1+q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2])-2q1 q3 Cos[t1]-2q2 q3 Cos[t2]+q3^2]]<= 10^-2]*1/Max[Abs[q3(W+Omega[q1]+Omega[q2])-(q1 Cos[t1]+ q2 Cos[t2])Omega[q3]],10^-100]*(1+1/Max[Abs[Exp[\[Beta] (W+Omega[q1]+Omega[q2]-Omega[q3])]-1],10^-100])*(1/(Exp[\[Beta] Omega[q1]]-1))(1/(Exp[\[Beta] Omega[q2]]-1))(1+1/(Exp[\[Beta] Omega[q3]]-1))*(q1^2 q2^2 q3^2 Sin[t1]Sin[t2])/(2Omega[q1]2Omega[q2])*1/4*(2Pi^2)/(2Pi)^9
myFullIntegrand[W_,\[Beta]_,q1_,q2_,t1_,t2_,phi_]:=Total[Map[myIntegrand[W,\[Beta],q1,q2,#,t1,t2,phi]&,DeltaRoots[W,\[Beta],q1,q2,t1,t2,phi]]]
myFullIntegrandCompiled=Compile[{W,\[Beta],q1,q2,t1,t2,phi},myFullIntegrand[W,\[Beta],q1,q2,t1,t2,phi]];
angleSubspaceInterpolation[W_,\[Beta]_,q1_,q2_]:=Interpolation[Flatten[Table[{{t1,t2,phi},myFullIntegrandCompiled[W,\[Beta],q1,q2,t1,t2,phi]},{t1,0,Pi,Pi/12},{t2,0,Pi,Pi/12},{phi,0,Pi,Pi/12}],2],InterpolationOrder->1]
angleSubspaceIntegral[W_,\[Beta]_,q1_,q2_]:=2*NIntegrate[angleSubspaceInterpolation[W,\[Beta],q1,q2][t1,t2,phi],{t1,0,Pi},{t2,0,Pi},{phi,0,Pi},PrecisionGoal->6, AccuracyGoal->30]
angleSubspaceIntegralCompiled=Compile[{W,\[Beta],q1,q2},angleSubspaceIntegral[W,\[Beta],q1,q2]];
momentumSubspaceInterpolation[W_,\[Beta]_]:=Interpolation[DeleteDuplicates[Flatten[Table[{{{q1,q2},#},{{q2,q1},#}}&[angleSubspaceIntegralCompiled[W,\[Beta],q1,q2]],{q1,0,Abs[W]+6,(Abs[W]+6)/50},{q2,0,q1,(Abs[W]+6)/50}],2]],InterpolationOrder->1]
workPoint[W_,\[Beta]_]:=2*NIntegrate[momentumSubspaceInterpolation[W,\[Beta]][q1,q2],{q1,0,Abs[W]+6},{q2,0,q1},PrecisionGoal->6,AccuracyGoal->30]
workPointCompiled = Compile[{W,\[Beta]},workPoint[W,\[Beta]]];
Print["End Initialization: "  DateString[Now]]
LaunchKernels[7];
distributionB17 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-.401,-.301,-.201,-.101,-.001,.001,.101}}];
CloseKernels[];
Export["/path/NtoNdist-B17.m",distributionB17];
Print["Batch #17: "  DateString[Now]]
LaunchKernels[7];
distributionB18 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{.201,.301,.401,.501,.601,.701,.801}}];
CloseKernels[];
Export["/path/NtoNdist-B18.m",distributionB18];
Print["Batch #18: "  DateString[Now]]
LaunchKernels[7];
distributionB19 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{.901,1.001,1.101,1.201,1.301,1.401,1.501}}];
CloseKernels[];
Export["/path/NtoNdist-B19.m",distributionB19];
Print["Batch #19: "  DateString[Now]]
LaunchKernels[7];
distributionB20 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{1.601,1.701,1.801,1.901,2.001,2.101,2.201}}];
CloseKernels[];
Export["/path/NtoNdist-B20.m",distributionB20];
Print["Batch #20: "  DateString[Now]]
LaunchKernels[7];
distributionB21 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{2.301,2.401,2.501,2.601,2.701,2.801,2.901}}];
CloseKernels[];
Export["/path/NtoNdist-B21.m",distributionB21];
Print["Batch #21: "  DateString[Now]]
LaunchKernels[7];
distributionB22 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{3.001,3.101,3.201,3.301,3.401,3.501,3.601}}];
CloseKernels[];
Export["/path/NtoNdist-B22.m",distributionB22];
Print["Batch #22: "  DateString[Now]]
LaunchKernels[7];
distributionB23 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{3.701,3.801,3.901,4.001,4.101,4.201,4.301}}];
CloseKernels[];
Export["/path/NtoNdist-B23.m",distributionB23];
Print["Batch #23: "  DateString[Now]]
LaunchKernels[7];
distributionB24 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{4.401,4.501,4.601,4.701,4.801,4.901,5.001}}];
CloseKernels[];
Export["/path/NtoNdist-B24.m",distributionB24];
Print["Batch #24: "  DateString[Now]]
LaunchKernels[7];
distributionB25 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{5.201,5.401,5.601,5.801,6.001,6.201,6.401}}];
CloseKernels[];
Export["/path/NtoNdist-B25.m",distributionB25];
Print["Batch #25: "  DateString[Now]]
LaunchKernels[7];
distributionB26 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{6.601,6.801,7.001,7.334,7.667,8.001,8.334}}];
CloseKernels[];
Export["/path/NtoNdist-B26.m",distributionB26];
Print["Batch #26: "  DateString[Now]]
LaunchKernels[7];
distributionB27 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{8.667,9.001,9.334,9.667,10.001,10.501,11.001}}];
CloseKernels[];
Export["/path/NtoNdist-B27.m",distributionB27];
Print["Batch #27: "  DateString[Now]]
LaunchKernels[7];
distributionB28 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{11.501,12.001,12.501,13.001,13.501,14.001,14.501}}];
CloseKernels[];
Export["/path/NtoNdist-B28.m",distributionB28];
Print["Batch #28: "  DateString[Now]]
LaunchKernels[7];
distributionB29 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{15.001,15.501,16.001,16.501,17.001,17.501,18.001}}];
CloseKernels[];
Export["/path/NtoNdist-B29.m",distributionB29];
Print["Batch #29: "  DateString[Now]]
LaunchKernels[7];
distributionB30 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{18.501,19.001,19.501,20.001,20.501,21.001,21.501}}];
CloseKernels[];
Export["/path/NtoNdist-B30.m",distributionB30];
Print["Batch #30: "  DateString[Now]]
LaunchKernels[7];
distributionB31 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{22.001,22.501,23.001,23.501,24.001,24.501,25.001}}];
CloseKernels[];
Export["/path/NtoNdist-B31.m",distributionB31];
Print["Batch #31: "  DateString[Now]]
LaunchKernels[7];
distributionB32 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{25.501,26.001,26.501,27.001,27.501,28.001,28.501}}];
CloseKernels[];
Export["/path/NtoNdist-B32.m",distributionB32];
Print["Batch #32: "  DateString[Now]]
LaunchKernels[3];
distributionB33 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{29.001,29.501,30.001}}];
CloseKernels[];
Export["/path/NtoNdist-B33.m",distributionB33];
Print["Batch #33: "  DateString[Now]]
Print["End Time: "  DateString[Now]]
