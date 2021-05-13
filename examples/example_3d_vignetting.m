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
hiaccuracy = 8; % for the large filter

% Light properties
wavelengths=linspace(0.73,0.85,100); % µm
polarization = 'unpolarized';

% Lens properties
fnumbers=[1.4];
cradegs=[0 5 10 20]; % chief ray angles


%% Create dielectric Fabry Perot filter using two materials


width=5.5; %micron

%% Define equivlent monolayer filter
% Target central wavelength
targetcwl = 0.800; %micron

% Incident medium and substrate
nair=1;
nsub=3.56; %silicon substarte

% Normalized bandwidth
normalized_fwhm=0.0130;

% Effective refractive index
effective_index=1.6;

% Equivalent filter
filter=tinyfilterCreateEquivalent(targetcwl,normalized_fwhm,effective_index,width,nair,nsub);
filterlarge=filter;         
filterlarge.width=100;
%% Pixel
pixel = pixel3D('width',width);
largepixel = pixel3D('width',100);


%% Run simulation for each fnumber and chief ray angle
for f=1:numel(fnumbers)
    for a=1:numel(cradegs)
        cradeg=cradegs(a);
        conedeg = atand(1./(2*fnumbers(f)));

        %% Simulate
        disp(['Simulate tiny filter: f/' num2str(fnumbers(f)) ' - CRA = ' num2str(cradeg) ' deg']);
        
        % Run simulation focusde light for equivalent monolayer
        azimuth_deg=0;
        eo16=load('lens-eo16.mat').lens;
        
        incident_wavepacket =  wavepacket3DLens(conedeg,cradeg,azimuth_deg);
        incident_wavepacket =  wavepacket3DLensVignetted(conedeg,cradeg,azimuth_deg,eo16.P,eo16.h,eo16.exitpupil);      
        
        Ttiny(:,a,f)=transmittanceTiny3D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixel);
        Tray(:,a,f) = transmittanceTinyRayFocused(lensVignetted(eo16.exitpupil,eo16.P,eo16.h),nair,effective_index,nsub,fwhm2reflectance(normalized_fwhm),width,targetcwl,wavelengths,conedeg,cradeg,polarization,accuracy,true);


       wavepacketIdealLens =  wavepacket3DLens(conedeg,cradeg,azimuth_deg);
        Tinf(:,a,f)=transmittanceTiny3D(filter,wavepacketIdealLens,wavelengths,polarization,hiaccuracy,pixel);
        
                
    end
end



%% Plot transmittance
cmap = hot;
s=size(cmap,1);
color{1}=cmap(1,:);
color{2}=cmap(round(0.4*s),:);
color{3}=cmap(round(0.6*s),:)
color{4}=cmap(round(0.65*s),:)

maxnorm = @(x) x;
count=1;

  fig=figure(f);clf;  hold on;
    fig.Position= [533 488 666 209];
for f =1:numel(fnumbers)
    
  

    for a=1:numel(cradegs)
        subplot(numel(fnumbers),numel(cradegs),count); hold on;
        hclassic(a)=plot(wavelengths,maxnorm(Tinf(:,a,f)),':','color','k','linewidth',1)
        htiny(a)=plot(wavelengths,maxnorm(Ttiny(:,a,f)),'color','r','linewidth',1)
        hray(a)=plot(wavelengths,maxnorm(Tray(:,a,f)),'color','b','linewidth',1)
        
        
        ylabel('Transmittance')
        xlabel('Wavelength (µm)')
        xlim([0.72 0.82])
        ylim([0 1])
        count=count+1;
    end
    title(['Filter transmittance for f/' num2str(fnumbers(f))])

    legend([htiny(1) hray(1) hclassic(1)],'Tiny Wave','Tiny Ray','Infinite filter','location','best')
end

