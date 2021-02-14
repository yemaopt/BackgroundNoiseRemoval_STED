% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

clear
clc
close all

addpath ./functions
addpath ./Data

z=-2:1/80:2;
L=length(z);
PatternExc=zeros(L,L,L);
PatternXYSTEDLCP=zeros(L,L,L);
PatternXYSTEDRCP=zeros(L,L,L);
PatternZSTED=zeros(L,L,L);
PatternEmi=zeros(L,L,L);

%% Excitation spot with (left) circular polarization
Intensity.Name='Gaussian';
Intensity.BeamWidth=1/0.3066; %mm measured by using knife-edge method
PhasePlate.Name='Null';
Polarization='Left circular'; %has tested that left and right circular polarization lead to the same result
for index=1:length(z)
    z(index)
    [ I,X,Y ] = FocusedPatternXY_CZT(Intensity,PhasePlate,Polarization,z(index),2,1/80);
    PatternExc(:,:,index)=I;
end

%% XY-STED spot with left circular polarization
Intensity.Name='Gaussian';
Intensity.BeamWidth=1/0.3066; %mm
PhasePlate.Name='Vortex 0-2pi';
Polarization='Left circular';

z=-2:1/80:2;
for index=1:length(z)
    z(index)
    [ I,X,Y ] = FocusedPatternXY_CZT(Intensity,PhasePlate,Polarization,z(index),2,1/80);
    PatternXYSTEDLCP(:,:,index)=I;
end

%%  XY-STED spot with right circular polarization
Intensity.Name='Gaussian';
Intensity.BeamWidth=1/0.3066; %mm
PhasePlate.Name='Vortex 0-2pi';
Polarization='Right circular';

z=-2:1/80:2;
for index=1:length(z)
    z(index)
    [ I,X,Y ] = FocusedPatternXY_CZT(Intensity,PhasePlate,Polarization,z(index),2,1/80);
    PatternXYSTEDRCP(:,:,index)=I;
end

%% Z-STED spot with (left) circular polarization
Intensity.Name='Gaussian';
Intensity.BeamWidth=1/0.3066*100/60; %mm
PhasePlate.Name='Binary 0-pi';
PhasePlate.InnerRadius=2;
Polarization='Left circular'; %has tested that left and right circular polarization lead to the same result

z=-2:1/80:2;
for index=1:length(z)
    z(index)
    [ I,X,Y ] = FocusedPatternXY_CZT(Intensity,PhasePlate,Polarization,z(index),2,1/80);
    PatternZSTED(:,:,index)=I;
end

%% Emission PSF with random polarization
Intensity.Name='Secant';
PhasePlate.Name='Null';
Polarization='Left circular'; %has tested that Linear X + Linear Y equals left circular pattern

z=-2:1/80:2;
for index=1:length(z)
    disp(z(index));
    [ I,X,Y ] = FocusedPatternXY_CZT(Intensity,PhasePlate,Polarization,z(index),2,1/80);
    PatternEmi(:,:,index)=I;
end

%% Display the results
figure(101);subplot(121);imagesc(PatternExc(:,:,(L+1)/2));axis equal;axis off;subplot(122);imagesc(squeeze(PatternExc(:,(L+1)/2,:))');axis equal;axis off;
figure(102);subplot(121);imagesc(PatternXYSTEDLCP(:,:,(L+1)/2));axis equal;axis off;subplot(122);imagesc(squeeze(PatternXYSTEDLCP(:,(L+1)/2,:))');axis equal;axis off;
figure(103);subplot(121);imagesc(PatternXYSTEDRCP(:,:,(L+1)/2));axis equal;axis off;subplot(122);imagesc(squeeze(PatternXYSTEDRCP(:,(L+1)/2,:))');axis equal;axis off;
figure(104);subplot(121);imagesc(PatternZSTED(:,:,(L+1)/2));axis equal;axis off;subplot(122);imagesc(squeeze(PatternZSTED(:,(L+1)/2,:))');axis equal;axis off;
figure(105);subplot(121);imagesc(PatternEmi(:,:,(L+1)/2));axis equal;axis off;subplot(122);imagesc(squeeze(PatternEmi(:,(L+1)/2,:))');axis equal;axis off;
%% save the results
Z=X;
save('FocusedPattern.mat','PatternExc','PatternEmi','PatternZSTED','PatternXYSTEDLCP','PatternXYSTEDRCP','X','Y','Z');
