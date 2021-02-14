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
PowerXYSTED=40; %mW
PowerZSTED=20; %mW

IntensityExc_XZ_Vectorized=PowerExc.*IntensityExc1uW_XZ(:)';
Intensity3DSTEDLCP_XZ=PowerXYSTED.*IntensityXYSTEDLCP1mW_XZ+PowerZSTED.*IntensityZSTED1mW_XZ;
Intensity3DSTEDLCP_XZ_Vectorized=Intensity3DSTEDLCP_XZ(:)';
Intensity3DSTEDRCP_XZ=PowerXYSTED.*IntensityXYSTEDRCP1mW_XZ+PowerZSTED.*IntensityZSTED1mW_XZ;
Intensity3DSTEDRCP_XZ_Vectorized=Intensity3DSTEDRCP_XZ(:)';

[L,~]=size(IntensityExc1uW_XZ);

tic
PSFexc_3DSTEDLCP_XZ=DetectedFluorescence(IntensityExc_XZ_Vectorized,Intensity3DSTEDLCP_XZ_Vectorized);
PSFexc_3DSTEDLCP_XZ=reshape(PSFexc_3DSTEDLCP_XZ,L,L);
toc
tic
PSFexc_3DSTEDRCP_XZ=DetectedFluorescence(IntensityExc_XZ_Vectorized,Intensity3DSTEDRCP_XZ_Vectorized);
PSFexc_3DSTEDRCP_XZ=reshape(PSFexc_3DSTEDRCP_XZ,L,L);
toc
tic
PSFexc_3DSTEDOnly_XZ=DetectedFluorescence(zeros(size(Intensity3DSTEDLCP_XZ_Vectorized)),Intensity3DSTEDLCP_XZ_Vectorized);
PSFexc_3DSTEDOnly_XZ=reshape(PSFexc_3DSTEDOnly_XZ,L,L);
toc
tic
PSFexc_3DSTEDGated_XZ=DetectedFluorescence_GatedDet(IntensityExc_XZ_Vectorized,Intensity3DSTEDLCP_XZ_Vectorized);
PSFexc_3DSTEDGated_XZ=reshape(PSFexc_3DSTEDGated_XZ,L,L);
toc

X=-2:1/80:2;
X_Exc=X;
Y_Exc=X;
Z_Exc=X;
[XX_Exc,ZZ_Exc]=meshgrid(X_Exc,Z_Exc);
[XXX_Exc,YYY_Exc,ZZZ_Exc]=meshgrid(X_Exc,Y_Exc,Z_Exc);
RRR_Exc=sqrt(XXX_Exc.^2+YYY_Exc.^2);
PSFexc_3DSTEDLCP_3D=zeros(L,L,L);
PSFexc_3DSTEDRCP_3D=zeros(L,L,L);
PSFexc_3DSTEDOnly_3D=zeros(L,L,L);
PSFexc_3DSTEDGated_3D=zeros(L,L,L);

%Rotate XZ cross sections to form 3D PSFs
for index=1:L
    disp(index);
    PSFexc_3DSTEDLCP_3D(index,:,:)=interp2(XX_Exc,ZZ_Exc,PSFexc_3DSTEDLCP_XZ',squeeze(RRR_Exc(index,:,:)),squeeze(ZZZ_Exc(index,:,:)),'linear',0);
    PSFexc_3DSTEDRCP_3D(index,:,:)=interp2(XX_Exc,ZZ_Exc,PSFexc_3DSTEDRCP_XZ',squeeze(RRR_Exc(index,:,:)),squeeze(ZZZ_Exc(index,:,:)),'linear',0);
    PSFexc_3DSTEDOnly_3D(index,:,:)=interp2(XX_Exc,ZZ_Exc,PSFexc_3DSTEDOnly_XZ',squeeze(RRR_Exc(index,:,:)),squeeze(ZZZ_Exc(index,:,:)),'linear',0);
    PSFexc_3DSTEDGated_3D(index,:,:)=interp2(XX_Exc,ZZ_Exc,PSFexc_3DSTEDGated_XZ',squeeze(RRR_Exc(index,:,:)),squeeze(ZZZ_Exc(index,:,:)),'linear',0);
end

PSF3DSTEDLCP=PSFexc_3DSTEDLCP_3D.*PSFdet1AU_3D;
PSF3DSTEDRCP=PSFexc_3DSTEDRCP_3D.*PSFdet1AU_3D;
PSF3DSTEDOnly=PSFexc_3DSTEDOnly_3D.*PSFdet1AU_3D;
PSF3DSTEDGated=PSFexc_3DSTEDGated_3D.*PSFdet1AU_3D;

%Sample is a single layer bead sample represented by a 321*321*321 matrix,
%1250 beads in total
%Pixel size: 640nm/80=8nm

load Sample1250Beads.mat

fSample=fftshift(fftn(ifftshift(Sample)));
fPSF3DSTEDLCP=fftshift(fftn(ifftshift(PSF3DSTEDLCP)));
fPSF3DSTEDRCP=fftshift(fftn(ifftshift(PSF3DSTEDRCP)));
fPSF3DSTEDOnly=fftshift(fftn(ifftshift(PSF3DSTEDOnly)));
fPSF3DSTEDGated=fftshift(fftn(ifftshift(PSF3DSTEDGated)));

Img3DSTEDLCP=fftshift(ifftn(ifftshift(fSample.*fPSF3DSTEDLCP)));
Img3DSTEDRCP=fftshift(ifftn(ifftshift(fSample.*fPSF3DSTEDRCP)));
Img3DSTEDOnly=fftshift(ifftn(ifftshift(fSample.*fPSF3DSTEDOnly)));
Img3DSTEDGated=fftshift(ifftn(ifftshift(fSample.*fPSF3DSTEDGated)));

ImgSTED=Img3DSTEDLCP;
ImgsubSTED=Img3DSTEDLCP-Img3DSTEDOnly;
ImggSTED=Img3DSTEDGated;
ImgpsSTED=Img3DSTEDLCP-Img3DSTEDRCP;

ImgSTED_XY=squeeze(ImgSTED(:,1:(L+1)/2,(L+1)/2));
ImgSTED_XZ=squeeze(ImgSTED(164,1:(L+1)/2,:))';
ImgSTED_XY=ImgSTED_XY./max(max(ImgSTED_XY));
ImgSTED_XZ=ImgSTED_XZ./max(max(ImgSTED_XZ));
ImgSTED_Combined=zeros(L,L);
ImgSTED_Combined(:,1:(L+1)/2)=ImgSTED_XY;
ImgSTED_Combined(:,(L+1)/2:L)=ImgSTED_XZ(:,1:161);

ImgsubSTED_XY=squeeze(ImgsubSTED(:,1:(L+1)/2,(L+1)/2));
ImgsubSTED_XZ=squeeze(ImgsubSTED(164,1:(L+1)/2,:))';
ImgsubSTED_XY=ImgsubSTED_XY./max(max(ImgsubSTED_XY));
ImgsubSTED_XZ=ImgsubSTED_XZ./max(max(ImgsubSTED_XZ));
ImgsubSTED_Combined=zeros(L,L);
ImgsubSTED_Combined(:,1:(L+1)/2)=ImgsubSTED_XY;
ImgsubSTED_Combined(:,(L+1)/2:L)=ImgsubSTED_XZ(:,1:161);

ImggSTED_XY=squeeze(ImggSTED(:,1:(L+1)/2,(L+1)/2));
ImggSTED_XZ=squeeze(ImggSTED(164,1:(L+1)/2,:))';
ImggSTED_XY=ImggSTED_XY./max(max(ImggSTED_XY));
ImggSTED_XZ=ImggSTED_XZ./max(max(ImggSTED_XZ));
ImggSTED_Combined=zeros(L,L);
ImggSTED_Combined(:,1:(L+1)/2)=ImggSTED_XY;
ImggSTED_Combined(:,(L+1)/2:L)=ImggSTED_XZ(:,1:161);

ImgpsSTED_XY=squeeze(ImgpsSTED(:,1:(L+1)/2,(L+1)/2));
ImgpsSTED_XZ=squeeze(ImgpsSTED(164,1:(L+1)/2,:))';
ImgpsSTED_XY=ImgpsSTED_XY./max(max(ImgpsSTED_XY));
ImgpsSTED_XZ=ImgpsSTED_XZ./max(max(ImgpsSTED_XZ));
ImgpsSTED_Combined=zeros(L,L);
ImgpsSTED_Combined(:,1:(L+1)/2)=ImgpsSTED_XY;
ImgpsSTED_Combined(:,(L+1)/2:L)=ImgpsSTED_XZ(:,1:161);

figure(101);imagesc(ImgSTED_Combined);colormap(hot(256));caxis([0,1]);axis equal;axis off;title('Regular STED')
figure(102);imagesc(ImgsubSTED_Combined);colormap(hot(256));caxis([0,1]);axis equal;axis off;title('sub-STED')
figure(103);imagesc(ImggSTED_Combined);colormap(hot(256));caxis([0,1]);axis equal;axis off;title('g-STED')
figure(104);imagesc(ImgpsSTED_Combined);colormap(hot(256));caxis([0,1]);axis equal;axis off;title('psSTED')

