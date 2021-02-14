% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

function y = Is(t,Isp,taus,deltat)

y=Isp.*exp(-(t-2.*taus-deltat).^2/(taus./2).^2);

end

