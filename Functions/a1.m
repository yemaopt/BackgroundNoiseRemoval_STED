% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

function [ y_a1 ] = a1( t, parameters )
c=299792458; %m/s
h=6.62606896e-34;

ksted=Is(t,parameters.Isp,parameters.taus,parameters.deltat).*parameters.sigma_depletion./(h*c/parameters.lambdas);
ksteddexc=Is(t,parameters.Isp,parameters.taus,parameters.deltat).*parameters.sigma_steddexc./(h*c/parameters.lambdas);
kexc=Ie(t,parameters.Iep,parameters.taue,parameters.taus).*parameters.sigma_exc./(h*c/parameters.lambdae);

kvib=1/parameters.tauv;
ks1=1/parameters.tauf;

gamma=ksted.*kvib./ks1./(kvib+ksted);
k1=ks1.*(1+gamma);
k2=ksteddexc+kexc;
k3=k2.*ksted./(kvib+ksted);

y_a1=-(k1+k2+k3);

end
