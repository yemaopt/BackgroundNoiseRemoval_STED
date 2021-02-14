% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

clear
clc
close all
addpath ./Functions
addpath ./Data

load FocusedPattern.mat

%Patterns are caculated as 321*321*321 matrices, with normalized coordinate
%of -2lambda:lambda/80:2lambda along each dimension

lambdaExc=640e-7; %cm
lambdaSTED=780e-7; %cm
T=1/(80e6); %Laser period, unit: second
TauExc=100e-12; %s
TauSTED=300e-12; %s
PowerExc1uW=1e-6;  %W
PowerSTED1mW=1e-3; %W
StepSize=lambdaExc/80; 

X_Exc=X.*lambdaExc;
Y_Exc=X_Exc;
Z_Exc=X_Exc;
X_STED=X.*lambdaSTED;
Y_STED=X_STED;
Z_STED=X_STED;

[XX_Exc,YY_Exc,ZZ_Exc]=meshgrid(X_Exc,Y_Exc,Z_Exc);
[XX_STED,YY_STED,ZZ_STED]=meshgrid(X_STED,Y_STED,Z_STED);

%Rescale the pattern to -2lambdaExc:lambdaExc/80:2lambdaExc
%PatternExc=PatternExc;
PatternXYSTEDLCP=interp3(XX_STED,YY_STED,ZZ_STED,PatternXYSTEDLCP,XX_Exc,YY_Exc,ZZ_Exc,'spline');
PatternXYSTEDRCP=interp3(XX_STED,YY_STED,ZZ_STED,PatternXYSTEDRCP,XX_Exc,YY_Exc,ZZ_Exc,'spline');
PatternZSTED=interp3(XX_STED,YY_STED,ZZ_STED,PatternZSTED,XX_Exc,YY_Exc,ZZ_Exc,'spline');

IntensityExc1uW_XZ=Power2Intensity_XZ(PowerExc1uW,PatternExc,StepSize,T,TauExc);
IntensityXYSTEDLCP1mW_XZ=Power2Intensity_XZ(PowerSTED1mW,PatternXYSTEDLCP,StepSize,T,TauSTED);
IntensityXYSTEDRCP1mW_XZ=Power2Intensity_XZ(PowerSTED1mW,PatternXYSTEDRCP,StepSize,T,TauSTED);
IntensityZSTED1mW_XZ=Power2Intensity_XZ(PowerSTED1mW,PatternZSTED,StepSize,T,TauSTED);

save Exc&STEDIntensity_XZ.mat...
     IntensityExc1uW_XZ...
     IntensityXYSTEDRCP1mW_XZ...
     IntensityXYSTEDLCP1mW_XZ...
     IntensityZSTED1mW_XZ

