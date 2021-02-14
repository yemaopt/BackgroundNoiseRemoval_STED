% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

clc
clear

addpath ./Functions
addpath ./Data
load Exc&STEDIntensity_XZ.mat
load PSFdet_1AU_3D.mat

PowerExc=4; %uW
IntensityExc_XZ_Vectorized=PowerExc.*IntensityExc1uW_XZ(:)';

[L,~]=size(IntensityExc1uW_XZ);

tic
PSFexc_Confocal_XZ=DetectedFluorescence(IntensityExc_XZ_Vectorized,zeros(size(IntensityExc_XZ_Vectorized)));
PSFexc_Confocal_XZ=reshape(PSFexc_Confocal_XZ,L,L);
toc

X=-2:1/80:2;
X_Exc=X;
Y_Exc=X;
Z_Exc=X;
[XX_Exc,ZZ_Exc]=meshgrid(X_Exc,Z_Exc);
[XXX_Exc,YYY_Exc,ZZZ_Exc]=meshgrid(X_Exc,Y_Exc,Z_Exc);
RRR_Exc=sqrt(XXX_Exc.^2+YYY_Exc.^2);
PSFexc_Confocal_3D=zeros(L,L,L);

%Rotate XZ cross section to form 3D PSF
for index=1:L
    disp(index);
    PSFexc_Confocal_3D(index,:,:)=interp2(XX_Exc,ZZ_Exc,PSFexc_Confocal_XZ',squeeze(RRR_Exc(index,:,:)),squeeze(ZZZ_Exc(index,:,:)),'linear',0);
end

PSFConfocal=PSFexc_Confocal_3D.*PSFdet1AU_3D;

%Sample is a single layer bead sample represented by a 321*321*321 matrix,
%1250 beads in total

load Sample1250Beads.mat

fSample=fftshift(fftn(ifftshift(Sample)));
fPSFConfocal=fftshift(fftn(ifftshift(PSFConfocal)));
ImgConfocal=fftshift(ifftn(ifftshift(fSample.*fPSFConfocal)));

ImgConfocal_XY=squeeze(ImgConfocal(:,1:(L+1)/2,(L+1)/2));
ImgConfocal_XZ=squeeze(ImgConfocal(164,1:(L+1)/2,:))';
ImgConfocal_XY=ImgConfocal_XY./max(max(ImgConfocal_XY));
ImgConfocal_XZ=ImgConfocal_XZ./max(max(ImgConfocal_XZ));
ImgConfocal_Combined=zeros(L,L);
ImgConfocal_Combined(:,1:(L+1)/2)=ImgConfocal_XY;
ImgConfocal_Combined(:,(L+1)/2:L)=ImgConfocal_XZ(:,1:161);

figure;imagesc(ImgConfocal_Combined);colormap(hot(256));caxis([0,1]);axis equal;axis off;title('Confocal')

