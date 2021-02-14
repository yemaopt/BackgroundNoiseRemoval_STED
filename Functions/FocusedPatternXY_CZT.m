% This code calculate the focused intensity pattern at the plane z (z=0 is
% the focal plane) using the vector diffraction theory
% The calculation is based on the chirp-z transform and the details are
% presented in the following reference:
% Marcel Leutenegger, Ramachandra Rao, Rainer A. Leitgeb, and Theo Lasser, 
% "Fast focus field calculations," Opt. Express 14, 11277-11291 (2006)

% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

function [ I,X,Y ] = FocusedPatternXY_CZT(Intensity,PhasePlate,Polarization,z,CalculationRange2,CalculationStep2)

lambda=1;
z2=z*lambda;
SamplingRange2=80*lambda;
SamplingStep2=1/80*lambda;
CalculationRange2=CalculationRange2*lambda;
CalculationStep2=CalculationStep2*lambda;

NA=1.4;
n=1.515;
theta3_max=asin(NA/n);
f=2; %mm The Focal length of the objective
PupilRadius=f*n*sin(theta3_max);
k0=2*pi/lambda;

Nfft=2*SamplingRange2/SamplingStep2+1;
Delta_kx=2*pi/SamplingStep2/Nfft;
Delta_ky=Delta_kx;

kx=-PupilRadius*k0/f:Delta_kx:PupilRadius*k0/f;
ky=-PupilRadius*k0/f:Delta_ky:PupilRadius*k0/f;
L=length(kx);
[kkx,kky]=meshgrid(kx,ky);

xx1=kkx./k0.*f;
yy1=kky./k0.*f;
phi3=atan2(yy1,xx1); %phi3=phi1
theta3=asin(sqrt(xx1.^2+yy1.^2)./f./n); %r1=f*n*sin(theta3)
Pupil=double((xx1.^2+yy1.^2)<=PupilRadius.^2);
theta3=theta3.*Pupil;

switch Intensity.Name
    case 'Uniform'
        Amp1=1.*Pupil;
    case 'Gaussian'
        w=Intensity.BeamWidth;
        Amp1=exp(-(xx1.^2+yy1.^2)./w^2).*Pupil;
    case 'Secant'
        Amp1=1./cos(theta3).*Pupil;
    otherwise
        error('No intensity definition');
end

switch PhasePlate.Name
    case 'Null'
        Phase=1;
    case 'Vortex 0-2pi'
        Phase=exp(1i.*phi3);
    case 'Binary 0-pi'
        rc=PhasePlate.InnerRadius;
        thetac=asin(rc/f/n);
        Phase=exp(1i.*pi.*double(theta3<=thetac));
        Phase=Phase.*Pupil;
    otherwise
        error('No phase plate definition');
end

switch Polarization
    case 'Left circular'
        E1=[1,-1i,0]'./sqrt(2);
    case 'Right circular'
        E1=[1,1i,0]'./sqrt(2);
    case 'Linear X'
        E1=[1,0,0]';
    case 'Linear Y'
        E1=[0,1,0]';
    otherwise
        error('No polarization definition');
end

%The following calculation is based on the vector diffraction theory
phi3=phi3(:);
theta3=theta3(:);
V1=[1+(cos(theta3)-1).*(cos(phi3).^2),(cos(theta3)-1).*cos(phi3).*sin(phi3),-sin(theta3).*cos(phi3)];
V2=[(cos(theta3)-1).*cos(phi3).*sin(phi3),1+(cos(theta3)-1).*(sin(phi3).^2),-sin(theta3).*sin(phi3)];
V3=[sin(theta3).*cos(phi3),sin(theta3).*sin(phi3),cos(theta3)];
E1x_prime=reshape(V1*E1,L,L);
E1y_prime=reshape(V2*E1,L,L);
E1z_prime=reshape(V3*E1,L,L);
theta3=reshape(theta3,L,L);

f1=1./sqrt(cos(theta3)).*E1x_prime.*exp(1i.*n.*k0.*z2.*cos(theta3)).*Amp1.*Phase./(f*n)^2; %sqrt(cos(theta3))/cos(theta3)
f2=1./sqrt(cos(theta3)).*E1y_prime.*exp(1i.*n.*k0.*z2.*cos(theta3)).*Amp1.*Phase./(f*n)^2;
f3=1./sqrt(cos(theta3)).*E1z_prime.*exp(1i.*n.*k0.*z2.*cos(theta3)).*Amp1.*Phase./(f*n)^2;

fs = 1/(Delta_kx/2/pi); freq1 = -CalculationRange2; freq2_0 = CalculationRange2 ;  
m = (freq2_0-freq1)/CalculationStep2+1;
freq2=m/(m-1)*(freq2_0-freq1)+freq1;
w = exp(-1j*2*pi*(freq2-freq1)/(m*fs)); 
a = exp(1j*2*pi*freq1/fs); 
    
E2x=Delta_kx^2.*czt(czt(f1,m,w,a)',m,w,a)';
E2y=Delta_kx^2.*czt(czt(f2,m,w,a)',m,w,a)';
E2z=Delta_kx^2.*czt(czt(f3,m,w,a)',m,w,a)';

I2x=abs(E2x).^2;
I2y=abs(E2y).^2;
I2z=abs(E2z).^2;

I=I2x+I2y+I2z;

X=((0:m-1)'*(freq2-freq1)/m) + freq1;
Y=X;

end

