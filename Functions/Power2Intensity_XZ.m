% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab,
% Department of Biomedical Engineering, Johns Hopkins University
% Last modified: Oct 25th, 2017

function [ IntensityXZ ] = Power2Intensity_XZ( Power,Pattern,step,T,tau )

Pattern=Pattern./max(max(max(Pattern)));
[N,~,~]=size(Pattern);
Scaler=1/(sum(sum(Pattern(:,:,(N+1)/2))).*step.*step);
Intensity3D=Pattern.*Power.*Scaler.*T./tau.*(2/sqrt(pi));
IntensityXZ=squeeze(Intensity3D(:,(N+1)/2,:));

end

