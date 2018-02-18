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
distributionB01 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-30.001,-29.501,-29.001,-28.501,-28.001,-27.501,-27.001}}];
CloseKernels[];
Export["/path/NtoNdist-B01.m",distributionB01];
Print["Batch #01: "  DateString[Now]]
LaunchKernels[7];
distributionB02 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-26.501,-26.001,-25.501,-25.001,-24.501,-24.001,-23.501}}];
CloseKernels[];
Export["/path/NtoNdist-B02.m",distributionB02];
Print["Batch #02: "  DateString[Now]]
LaunchKernels[7];
distributionB03 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-23.001,-22.501,-22.001,-21.501,-21.001,-20.501,-20.001}}];
CloseKernels[];
Export["/path/NtoNdist-B03.m",distributionB03];
Print["Batch #03: "  DateString[Now]]
LaunchKernels[7];
distributionB04 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-19.501,-19.001,-18.501,-18.001,-17.501,-17.001,-16.501}}];
CloseKernels[];
Export["/path/NtoNdist-B04.m",distributionB04];
Print["Batch #05: "  DateString[Now]]
LaunchKernels[7];
distributionB06 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-16.001,-15.501,-15.001,-14.501,-14.001,-13.501,-13.001}}];
CloseKernels[];
Export["/path/NtoNdist-B06.m",distributionB06];
Print["Batch #06: "  DateString[Now]]
LaunchKernels[7];
distributionB07 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-12.501,-12.001,-11.501,-11.001,-10.501,-10.001,-9.667}}];
CloseKernels[];
Export["/path/NtoNdist-B07.m",distributionB07];
Print["Batch #07: "  DateString[Now]]
LaunchKernels[7];
distributionB08 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-9.334,-9.001,-8.667,-8.334,-8.001,-7.667,-7.334}}];
CloseKernels[];
Export["/path/NtoNdist-B08.m",distributionB08];
Print["Batch #08: "  DateString[Now]]
LaunchKernels[7];
distributionB09 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-7.001,-6.801,-6.601,-6.401,-6.201,-6.001,-5.801}}];
CloseKernels[];
Export["/path/NtoNdist-B09.m",distributionB09];
Print["Batch #09: "  DateString[Now]]
LaunchKernels[7];
distributionB10 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-5.601,-5.401,-5.201,-5.001,-4.901,-4.801,-4.701}}];
CloseKernels[];
Export["/path/NtoNdist-B10.m",distributionB10];
Print["Batch #10: "  DateString[Now]]
LaunchKernels[7];
distributionB11 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-4.601,-4.501,-4.401,-4.301,-4.201,-4.101,-4.001}}];
CloseKernels[];
Export["/path/NtoNdist-B11.m",distributionB11];
Print["Batch #11: "  DateString[Now]]
LaunchKernels[7];
distributionB12 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-3.901,-3.801,-3.701,-3.601,-3.501,-3.401,-3.301}}];
CloseKernels[];
Export["/path/NtoNdist-B12.m",distributionB12];
Print["Batch #12: "  DateString[Now]]
LaunchKernels[7];
distributionB13 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-3.201,-3.101,-3.001,-2.901,-2.801,-2.701,-2.601}}];
CloseKernels[];
Export["/path/NtoNdist-B13.m",distributionB13];
Print["Batch #13: "  DateString[Now]]
LaunchKernels[7];
distributionB14 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-2.501,-2.401,-2.301,-2.001,-1.901,-1.801,-1.701}}];
CloseKernels[];
Export["/path/NtoNdist-B14.m",distributionB14];
Print["Batch #14: "  DateString[Now]]
LaunchKernels[7];
distributionB15 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-1.601,-1.501,-1.401,-2.201,-2.101,-1.301,-1.201}}];
CloseKernels[];
Export["/path/NtoNdist-B15.m",distributionB15];
Print["Batch #15: "  DateString[Now]]
LaunchKernels[7];
distributionB16 = ParallelTable[{W,workPointCompiled[W,temperature]},{W,{-1.101,-1.001,-.901,-.801,-.701,-.601,-.501}}];
CloseKernels[];
Export["/path/NtoNdist-B16.m",distributionB16];
Print["Batch #16: "  DateString[Now]]
Print["End Time: "  DateString[Now]]
