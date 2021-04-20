clear; close all;

addpath('/home/thomas/Documents/tinyfilters/research/wavepacket')
addpath(genpath('../../core'))

%% Create dielectric Fabry Perot filter using two materials

% Target central wavelength
targetcwl = 0.791; %micron


nair=1;
nsub=3.56; %silicon substarte

    
%% Choose simulation options

polarization = 'unpolarized';

accuracy = 8;
wavelengths=linspace(0.73,0.85,200); % µm




fnumber=8;
conedeg = atand(1./(2*fnumber));
cradegs=[0 5 10 20]


%% Choose filter properties
fwhm=9*1e-3; % µm

L =  fwhm/targetcwl;

width=5.5;
pwidth=sqrt(0.42)*width;
pixelkernel=pixel_fullwidth(width);
neff=1.7



%% Run simulation for each angle
for a=1:numel(cradegs)
    cradeg=cradegs(a);
    %% Simulate
    disp(['Simulate tiny filter: CRA = ' num2str(cradeg) ' deg']);
    
    
    Ttiny(:,a)=tinytransmittance3dfocus_rotational_mono(targetcwl,L,neff,width,nsub,conedeg,cradeg,0,wavelengths,polarization,pixelkernel,accuracy);
    
    Tclassic(:,a)=tinytransmittance3dfocus_rotational_mono(targetcwl,L,neff,100,nsub,conedeg,cradeg,0,wavelengths,polarization,pixel_fullwidth(100),accuracy);
    
end




%% Plot transmittance
% There there is a drop in transmittance and an increase in FWHM
% 

cmap = hot;
s=size(cmap,1);
color{1}=cmap(1,:);
color{2}=cmap(round(0.4*s),:);
color{3}=cmap(round(0.6*s),:)
color{4}=cmap(round(0.65*s),:)


addpath('/home/thomas/Documents/imec/phd/scripts/fabryperotconvolution/kernels')


fig=figure(1);clf;  hold on;
fig.Position= [533 488 666 209];
for a=1:numel(cradegs)
    hclassic(a)=plot(wavelengths,Tclassic(:,a),':','color',color{a},'linewidth',1)
    htiny(a)=plot(wavelengths,Ttiny(:,a),'color',color{a},'linewidth',1)
    
    ylabel('Transmittance')
    xlabel('Wavelength (µm)')
    xlim([0.72 0.82])
    ylim([0 1])
end
title(['Filter transmittance for f/' num2str(fnumber)])

legend([htiny hclassic(1)],'CRA = 0 deg','CRA = 5 deg','CRA = 10 deg','CRA = 20 deg','Infinite filter','location','best')