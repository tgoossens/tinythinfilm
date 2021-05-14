    clear;
    
 %% Comparison of tiny Fabry-Pérot filter and equivalent tiny monolayer
% This example demonstrated how the transmittance of a tiny dielectric Fabry-Pérot filter
% differs from equivalent monolayer response.
%
% It shows that a well chosen equivalent monolayer approximates the response of a Fabry-pérot well. That is, given that the Fabry-Pérot is nicely approximated by a Lorentzian.
%
%
% Copyright Thomas Goossens

clear; close all;


addpath('/home/thomas/Documents/tinyfilters/research/wavepacket')

%% Create dielectric Fabry Perot filter using two materials

% Target central wavelength (where the filter is centered)
targetcwl = 0.800; %micron

filterwidth=2;
%% Refractive indeices of indicent medium and substrate
nair=1;
nsub=3.56; %silicon substrate


%% Tiny Fabry-Pérot filter design
neff=1.7


%% Equivalent filter
normalized_fwhm=0.01;
filter=tinyfilterCreateEquivalent(targetcwl,normalized_fwhm,neff,filterwidth,nair,nsub);
pixel=pixel2D('range',[0 0.5*filterwidth])

wavepacket=wavepacket2DCollimated(10,nair)
displaySetup(filter,'pixel',pixel,'wavepacket',wavepacket)


%%
