%%
clear; 



%% Load data crosstlak
%data=load('data-crosstalk.mat');





%%
n0=1;
neff=1.71;
nsub=1;
R=[fwhm2reflectance(0.007)];
Rs=[R R];
filterwidths=[5*5.5 5.5];
cwls_micron = [713.8 720]/1000;
wavelengths_micron=(600:800)/1000;
polarization='s'
accuracy = 10;
pixelrange=0.5*[-filterwidths(2) filterwidths(2)]

largepixelrange=[0 filterwidths];
% 
% % Fit a at normal incidence
% angledeg=0.01;
% Tref_right=transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidths(2),cwls_micron(2),wavelengths_micron,angledeg,polarization,accuracy,pixelrange,true)
% Tneighbour=transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,2*filterwidths(1),cwls_micron(1),wavelengths_micron,angledeg,polarization,accuracy,largepixelrange,true)
% 
% figure;clf; hold on
% plot(data.w/1000,data.measured/max(data.measured),'k');
% plot(data.w/1000,data.neighbour/max(data.neighbour),'m--');
% plot(wavelengths_micron,Tref_right,'r')
% plot(wavelengths_micron,Tneighbour,'b')
% xlim([0.6 0.7])
%%
angles=[0.1 5 10 15 20]
fig=figure;clf;hold on;
fig.Position=[389 365 1521 407];
fig.Position=[1238 249 471 622];
for a=1:numel(angles)
subplot(numel(angles),1,a); hold on;    
angledeg=angles(a);

T=transmittanceTinyRayEquivalent_crosstalk(n0,neff,nsub,Rs,filterwidths,cwls_micron,wavelengths_micron,angledeg,polarization,accuracy,pixelrange)
Tref_right=transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidths(2),cwls_micron(2),wavelengths_micron,angledeg,polarization,accuracy,pixelrange,true)
Tref_left=transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,2*filterwidths(1),cwls_micron(1),wavelengths_micron,angledeg,polarization,accuracy,largepixelrange,true)


plot(wavelengths_micron,(T),'r')
plot(wavelengths_micron,Tref_right,'r:')
plot(wavelengths_micron,(Tref_left),'k')
xlim([600 750]/1000)
%ylim([0 1])
title(num2str(angledeg))
end