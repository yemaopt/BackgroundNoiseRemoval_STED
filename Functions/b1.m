% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

function [ y_b1 ] = b1( t, parameters )
c=299792458; %m/s
h=6.62606896e-34;

ksteddexc=Is(t,parameters.Isp,parameters.taus,parameters.deltat).*parameters.sigma_steddexc./(h*c/parameters.lambdas);
kexc=Ie(t,parameters.Iep,parameters.taue,parameters.taus).*parameters.sigma_exc./(h*c/parameters.lambdae);

k2=ksteddexc+kexc;

y_b1=k2;

end

