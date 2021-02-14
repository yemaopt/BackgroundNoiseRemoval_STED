% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

clc
clear
close all

addpath ./Functions

Isp=0:1e6:5e9;
Iep=6.6321e+05.*ones(1,length(Isp)); %W/cm2 -> 4uW

deltat=0; %s
F_0ps=DetectedFluorescence_OptDeltat(Iep,Isp,deltat);
deltat=50E-12; %s
F_50ps=DetectedFluorescence_OptDeltat(Iep,Isp,deltat);
deltat=150E-12; %s
F_150ps=DetectedFluorescence_OptDeltat(Iep,Isp,deltat);
deltat=250E-12; %s
F_250ps=DetectedFluorescence_OptDeltat(Iep,Isp,deltat);

figure(101);
hold on
plot(Isp,log10(F_0ps./max(F_0ps)))
plot(Isp,log10(F_50ps./max(F_50ps)))
plot(Isp,log10(F_150ps./max(F_150ps)))
plot(Isp,log10(F_250ps./max(F_250ps)))
xlabel('STED Intensity (W/cm^2)')
ylabel('Log_1_0(Normalized fluorescence)')
legend('0ps','50ps','150ps','250ps')

