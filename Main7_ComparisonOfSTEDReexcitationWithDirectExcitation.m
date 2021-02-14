% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

clc
clear
close all

addpath ./Functions

Isp=0:1e6:5e9;
Iep=6.6321e+05.*ones(1,length(Isp)); %W/cm2 -> 4uW

STEDDirExcFlag=1;
[F_wDirwRe,F_wDirwoRe]=DetectedFluorescence_TestSTEDReexcDirexc(Iep,Isp,STEDDirExcFlag);
STEDDirExcFlag=0;
[F_woDirwRe,F_woDirwoRe]=DetectedFluorescence_TestSTEDReexcDirexc(Iep,Isp,STEDDirExcFlag);

figure(102);
hold on
plot(Isp,log10(F_wDirwRe./max(F_wDirwRe)),'-.')
plot(Isp,log10(F_wDirwoRe./max(F_wDirwoRe)),':')
plot(Isp,log10(F_woDirwRe./max(F_woDirwRe)))
plot(Isp,log10(F_woDirwoRe./max(F_woDirwoRe)),'--')
xlabel('STED Intensity (W/cm^2)')
ylabel('Log_1_0(Normalized fluorescence)')
legend('w/ STED dir-exc, w/ STED re-exc','w/ STED dir-exc, w/o STED re-exc','w/o STED dir-exc, w/ STED re-exc','w/o STED dir-exc, w/o STED re-exc')




