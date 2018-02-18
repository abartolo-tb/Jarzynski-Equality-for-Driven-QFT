(* ::Package:: *)

Print["Start Time: "  DateString[Now]]
temperature=1;
Omega[q_]:=Sqrt[1+q^2]
myDiscriminant[W_,q1_,q2_,t1_,t2_,phi_]:= ((q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2])-(W-Omega[q1]-Omega[q2])^2)/Max[Abs[2(W-Omega[q1]-Omega[q2])],10^-100])^2+((q1 Cos[t1] + q2 Cos[t2])/Max[Abs[(W-Omega[q1]-Omega[q2])],10^-100])^2-1
DeltaRoots[W_,\[Beta]_,q1_,q2_,t1_,t2_,phi_]:=If[myDiscriminant[W,q1,q2,t1,t2,phi]<= 0,{10^-100},{1/(1-((q1 Cos[t1] + q2 Cos[t2])/Max[Abs[(W-Omega[q1]-Omega[q2])],10^-100])^2) (-(1/2)(1-(q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2]))/(Max[Abs[(W-Omega[q1]-Omega[q2])],10^-100])^2)(q1 Cos[t1] + q2 Cos[t2])+Sqrt[myDiscriminant[W,q1,q2,t1,t2,phi]]),1/(1-((q1 Cos[t1] + q2 Cos[t2])/Max[Abs[(W-Omega[q1]-Omega[q2])],10^-100])^2) (-(1/2)(1-(q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2]))/(Max[Abs[(W-Omega[q1]-Omega[q2])],10^-100])^2)(q1 Cos[t1] + q2 Cos[t2])-Sqrt[myDiscriminant[W,q1,q2,t1,t2,phi]])},{10^-100}]
myIntegrand[W_,\[Beta]_,q1_,q2_,q3_,t1_,t2_,phi_]:=Boole[q3>10^-50]*Boole[ Abs[W-Sqrt[1+q1^2]-Sqrt[1+q2^2]-Sqrt[1+q3^2]+Sqrt[1+q1^2+q2^2+2q1 q2 (Sin[t1]Sin[t2]Cos[phi]+Cos[t1]Cos[t2])+2q1 q3 Cos[t1]+2q2 q3 Cos[t2]+q3^2]]<= 10^-2]*1/Max[Abs[q3(W- Omega[q1]-Omega[q2])+(q1 Cos[t1]+ q2 Cos[t2])Omega[q3]],10^-100]*(1/Max[Abs[Exp[\[Beta] ( Omega[q1]+Omega[q2]+Omega[q3]-W)]-1],10^-100])*(1+1/(Exp[\[Beta] Omega[q1]]-1))(1+1/(Exp[\[Beta] Omega[q2]]-1))(1+1/(Exp[\[Beta] Omega[q3]]-1))*(q1^2 q2^2 q3^2 Sin[t1]Sin[t2])/(2Omega[q1]2Omega[q2])*1/3!*(2Pi^2)/(2Pi)^9
myFullIntegrand[W_,\[Beta]_,q1_,q2_,t1_,t2_,phi_]:=Total[Map[myIntegrand[W,\[Beta],q1,q2,#,t1,t2,phi]&,DeltaRoots[W,\[Beta],q1,q2,t1,t2,phi]]]
myFullIntegrandCompiled=Compile[{W,\[Beta],q1,q2,t1,t2,phi},myFullIntegrand[W,\[Beta],q1,q2,t1,t2,phi]];
angleSubspaceInterpolation[W_,\[Beta]_,q1_,q2_]:=Interpolation[Flatten[Table[{{t1,t2,phi},myFullIntegrandCompiled[W,\[Beta],q1,q2,t1,t2,phi]},{t1,0,Pi,Pi/12},{t2,0,Pi,Pi/12},{phi,0,Pi,Pi/12}],2],InterpolationOrder->1]
angleSubspaceIntegral[W_,\[Beta]_,q1_,q2_]:=2*NIntegrate[angleSubspaceInterpolation[W,\[Beta],q1,q2][t1,t2,phi],{t1,0,Pi},{t2,0,Pi},{phi,0,Pi},PrecisionGoal->6, AccuracyGoal->20]
angleSubspaceIntegralCompiled=Compile[{W,\[Beta],q1,q2},angleSubspaceIntegral[W,\[Beta],q1,q2]];
momentumSubspaceInterpolation[W_,\[Beta]_]:=Interpolation[DeleteDuplicates[Flatten[Table[{{{q1,q2},#},{{q2,q1},#}}&[angleSubspaceIntegralCompiled[W,\[Beta],q1,q2]],{q1,0,W+4,(W+4)/50},{q2,0,q1,(W+4)/50}],2]],InterpolationOrder->1]
workPoint[W_,\[Beta]_]:=2*NIntegrate[momentumSubspaceInterpolation[W,\[Beta]][q1,q2],{q1,0,W+4},{q2,0,q1},PrecisionGoal->6,AccuracyGoal->20]
workPointCompiled = Compile[{W,\[Beta]},workPoint[W,\[Beta]]];
Print["End Initialization: "  DateString[Now]]
LaunchKernels[7];
distributionB01 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{0.001,0.201,0.401,0.601,0.801,1.001,1.201}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B01.m",distributionB01];
Print["Batch #01: "  DateString[Now]]
LaunchKernels[7];
distributionB02 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{1.401,1.601,1.801,2.001,2.101,2.201,2.301}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B02.m",distributionB02];
Print["Batch #02: "  DateString[Now]]
LaunchKernels[7];
distributionB03 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{2.401,2.501,2.601,2.701,2.801,2.901,3.001}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B03.m",distributionB03];
Print["Batch #03: "  DateString[Now]]
LaunchKernels[7];
distributionB04 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{3.101,3.201,3.301,3.401,3.501,3.601,3.701}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B04.m",distributionB04];
Print["Batch #04: "  DateString[Now]]
LaunchKernels[7];
distributionB05 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{3.801,3.901,4.001,4.101,4.201,4.301,4.401}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B05.m",distributionB05];
Print["Batch #05: "  DateString[Now]]
LaunchKernels[7];
distributionB06 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{4.501,4.601,4.701,4.801,4.901,5.001,5.201}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B06.m",distributionB06];
Print["Batch #06: "  DateString[Now]]
LaunchKernels[7];
distributionB07 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{5.401,5.601,5.801,6.001,6.201,6.401,6.601}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B07.m",distributionB07];
Print["Batch #07: "  DateString[Now]]
LaunchKernels[7];
distributionB08 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{6.801,7.001,7.201,7.401,7.601,7.801,8.001}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B08.m",distributionB08];
Print["Batch #08: "  DateString[Now]]
LaunchKernels[7];
distributionB09 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{8.201,8.401,8.601,8.801,9.001,9.201,9.401}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B09.m",distributionB09];
Print["Batch #09: "  DateString[Now]]
LaunchKernels[7];
distributionB10 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{9.601,9.801,10.001,10.334,10.667,11.001,11.334}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B10.m",distributionB10];
Print["Batch #10: "  DateString[Now]]
LaunchKernels[7];
distributionB11 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{11.667,12.001,12.334,12.667,13.001,13.334,13.667}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B11.m",distributionB11];
Print["Batch #11: "  DateString[Now]]
LaunchKernels[7];
distributionB12 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{14.001,14.334,14.667,15.001,15.334,15.667,16.001}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B12.m",distributionB12];
Print["Batch #12: "  DateString[Now]]
LaunchKernels[7];
distributionB13 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{16.334,16.667,17.001,17.334,17.667,18.001,18.334}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B13.m",distributionB13];
Print["Batch #13: "  DateString[Now]]
LaunchKernels[7];
distributionB14 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{18.667,19.001,19.334,19.667,20.001,20.334,20.667}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B14.m",distributionB14];
Print["Batch #14: "  DateString[Now]]
LaunchKernels[7];
distributionB15 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{21.001,21.334,21.667,22.001,22.334,22.667,23.001}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B15.m",distributionB15];
Print["Batch #15: "  DateString[Now]]
LaunchKernels[7];
distributionB16 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{23.334,23.667,24.001,24.334,24.667,25.001,25.334}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B16.m",distributionB16];
Print["Batch #16: "  DateString[Now]]
LaunchKernels[7];
distributionB17 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{25.667,26.001,26.334,26.667,27.001,27.334, 27.667}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B17.m",distributionB17];
Print["Batch #17: "  DateString[Now]]
LaunchKernels[7];
distributionB18 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{28.001,28.334,28.667,29.001,29.334,29.667,30.001}}];
CloseKernels[];
Export["/path/NtoNplus2dist-B18.m",distributionB18];
Print["Batch #18: "  DateString[Now]]
Print["End Time: "  DateString[Now]]
