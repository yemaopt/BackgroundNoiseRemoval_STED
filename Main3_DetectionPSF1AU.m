% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

clc
clear
close all

addpath ./Functions
addpath ./Data

load FocusedPattern.mat

lambdaExc=640e-7; %cm
lambdaEmi=670e-7; %cm

X_Exc=X.*lambdaExc;
Y_Exc=X_Exc;
Z_Exc=X_Exc;
X_Emi=X.*lambdaEmi;
Y_Emi=X_Emi;
Z_Emi=X_Emi;

[XX_Exc,YY_Exc,ZZ_Exc]=meshgrid(X_Exc,Y_Exc,Z_Exc);
[XX_Emi,YY_Emi,ZZ_Emi]=meshgrid(X_Emi,Y_Emi,Z_Emi);

PatternEmi=interp3(XX_Emi,YY_Emi,ZZ_Emi,PatternEmi,XX_Exc,YY_Exc,ZZ_Exc,'spline');

[XX_Exc,YY_Exc]=meshgrid(X_Exc,Y_Exc);
PinholeDiameter=1.22*lambdaEmi/1.4;
Pinhole=double((XX_Exc.^2+YY_Exc.^2)<=(PinholeDiameter/2)^2);
fPinhole=fftshift(fft2(ifftshift(Pinhole)));

L=length(Z_Exc);
for i=1:L
    disp(i);
    fPatternEmi=fftshift(fft2(ifftshift(PatternEmi(:,:,i))));
    PSFdet1AU_3D(:,:,i)=fftshift(ifft2(ifftshift(fPatternEmi.*fPinhole)));
end

PSFdet1AU_3D=PSFdet1AU_3D./max(max(max(PSFdet1AU_3D)));
PSFdet0AU_3D=PatternEmi./max(max(max(PatternEmi)));

figure
imagesc(squeeze(PSFdet0AU_3D(:,(L+1)/2,:))')
axis off
axis equal

figure
imagesc(squeeze(PSFdet1AU_3D(:,(L+1)/2,:))')
axis off
axis equal

save PSFdet_1AU_3D.mat PSFdet1AU_3D


