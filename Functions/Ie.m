% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

function y = Ie(t,Iep,taue,taus)

y=Iep.*exp(-(t-2.*taus).^2/(taue./2).^2);

end

