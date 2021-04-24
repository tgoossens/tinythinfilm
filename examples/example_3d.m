%% Example of 3D pixel-integration Fabry-Pérot filter using an equivalent monolayer model
%
% This example demonstrates the effect of on-axis and off-axis illumination on the effective transmittance and shows the difference
% for an large vs tiny filter.
% 
%
% The simulation is performed for f/2.8 and f/8. 
% At f/8 the tiny-filter effect dominates the angular dependency while for f/2.8 the additional smoothing by the cone angle dominates.
% Copyright Thomas Goossens

clear; close all;


addpath('/home/thomas/Documents/tinyfilters/research/wavepacket')
addpath(genpath('../core'))

    
%% Choose simulation options

% Accuracy options
accuracy = 6;

% Light properties
wavelengths=linspace(0.73,0.85,100); % µm
polarization = 'unpolarized';

% Lens properties
fnumbers=[2.8 8]; 
cradegs=[0 5 10 20]; % chief ray angles


%% Create dielectric Fabry Perot filter using two materials

% Target central wavelength
targetcwl = 0.800; %micron


nair=1;
nsub=3.56; %silicon substarte

nl = 1.4; % low refractive index
nh = 2.4; % high refractive index

dh = targetcwl/(4*nh);%quarterwave 
dl = targetcwl/(4*nl);%quarterwave 

n = [nh nl nh nl nh nl nh [nl nl] nh nl nh nl nh nl nh];
thickness = [dh dl dh dl dh dl dh [dl dl] dh dl dh dl dh dl dh];

width=5.5; %micron


filter=tinyfilter(nair,n,nsub,thickness,width);
normalized_fwhm=0.0130;
effective_index=1.6;

pixelkernel = pixel_fullwidth(width);


%% Run simulation for each fnumber and chief ray angle
for f=1:numel(fnumbers)
    for a=1:numel(cradegs)
        cradeg=cradegs(a);
        conedeg = atand(1./(2*fnumbers(f)));

        %% Simulate
        disp(['Simulate tiny filter: f/' num2str(fnumbers(f)) ' - CRA = ' num2str(cradeg) ' deg']);
        
        % Run simulation focusde light for equivalent monolayer
        azimuth_deg=0;
        incident_wavepacket =  wavepacket3d_focus(conedeg,cradeg,azimuth_deg,width);
        
        Ttiny(:,a,f)=tinytransmittance3d_core(filter2,incident_wavepacket,wavelengths,polarization,accuracy,pixelkernel);
        
        large_wavepacket=  wavepacket3d_focus(conedeg,cradeg,azimuth_deg,100);
        Tclassic(:,a,f)=tinytransmittance3d_core(filter2,large_wavepacket,wavelengths,polarization,accuracy,pixel_fullwidth(100));
        
                
    end
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


for f =1:numel(fnumbers)
    fig=figure(f);clf;  hold on;
    fig.Position= [533 488 666 209];

    for a=1:numel(cradegs)
        hclassic(a)=plot(wavelengths,Tclassic(:,a,f),':','color',color{a},'linewidth',1)
        htiny(a)=plot(wavelengths,Ttiny(:,a,f),'color',color{a},'linewidth',1)
        htiny(a)=plot(wavelengths,Ttiny(:,a,f),'color',color{a},'linewidth',1)
        
        ylabel('Transmittance')
        xlabel('Wavelength (µm)')
        xlim([0.72 0.82])
        ylim([0 1])
    end
    title(['Filter transmittance for f/' num2str(fnumbers(f))])

    legend([htiny  hclassic(1)],'CRA = 0 deg','CRA = 5 deg','CRA = 10 deg','CRA = 20 deg','Infinite filter','location','best')
end

