% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

function [ Fluorescence, Fluorescence_woReexc ] = DetectedFluorescence_TestSTEDReexcDirexc( Iep,Isp,STEDDirExcFlag )

T=1.25e-8; %s
dt=2e-12; %s

parameters.sigma_exc=5.7e-16; %cm2
parameters.sigma_depletion=0.23e-16; %cm2
if STEDDirExcFlag==1
    parameters.sigma_steddexc=5.7e-21; %cm2
else
    parameters.sigma_steddexc=0; %cm2
end

parameters.tauf=3.5e-9; %s
parameters.tauv=0.2e-12; %s
parameters.lambdae=640e-9; %m
parameters.lambdas=780e-9; %m
parameters.taue=100e-12; %s
parameters.taus=300e-12; %s
parameters.deltat=150e-12; %s

parameters.Iep=Iep; %W/cm2 
parameters.Isp=Isp; %W/cm2 

t=dt:dt:T;
Ps1=zeros(length(t),length(parameters.Iep));
Ps1_woReexc=zeros(length(t),length(parameters.Iep));
Ps1_0=0;

for index=1:1:length(t)
    if index==1
        Ps1(index,:)=exp(dt.*a1(0,parameters)).*(Ps1_0+b1(0,parameters)./a1(0,parameters))-b1(0,parameters)./a1(0,parameters);
        Ps1_woReexc(index,:)=exp(dt.*a1_woReexc(0,parameters)).*(Ps1_0+b1(0,parameters)./a1_woReexc(0,parameters))-b1(0,parameters)./a1_woReexc(0,parameters);
    else
        Ps1(index,:)=exp(dt.*a1(t(index-1),parameters)).*(Ps1(index-1,:)+b1(t(index-1),parameters)./a1(t(index-1),parameters))-b1(t(index-1),parameters)./a1(t(index-1),parameters);
        Ps1_woReexc(index,:)=exp(dt.*a1_woReexc(t(index-1),parameters)).*(Ps1_woReexc(index-1,:)+b1(t(index-1),parameters)./a1_woReexc(t(index-1),parameters))-b1(t(index-1),parameters)./a1_woReexc(t(index-1),parameters);
    end
end
disp('Done');
SumPs1=sum(Ps1,1).*dt;
SumPs1_woReexc=sum(Ps1_woReexc,1).*dt;

QuantumYeild=0.65;
Tdwelling=400e-6; %s
yita=(2*pi*(1-cos(asin(1.4/1.515))))/(4*pi)*0.9; %APD effciency 90%

Fluorescence=SumPs1.*QuantumYeild.*yita.*Tdwelling./T./parameters.tauf;
Fluorescence_woReexc=SumPs1_woReexc.*QuantumYeild.*yita.*Tdwelling./T./parameters.tauf;
end

