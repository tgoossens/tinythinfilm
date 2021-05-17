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
accuracy = 7;
hiaccuracy = 8; % for the large filter

% Light properties
wavelengths=linspace(0.73,0.85,100); % µm
polarization = 'unpolarized';

% Lens properties
fnumbers=[2.8 1.4];
cradegs=[10 20]; % chief ray angles






%% Define equivlent monolayer filter

% Width
width=3.56 %micron

% Target central wavelength
targetcwl = 0.800; %micron

% Incident medium and substrate
nair=1;
nsub=3.56; %silicon substarte

% Normalized bandwidth
normalized_fwhm=0.0114;

% Effective refractive index
effective_index=1.7;

% Equivalent filter
filter=tinyfilterCreateEquivalent(targetcwl,normalized_fwhm,effective_index,width,nair,nsub);
filterlarge=filter;         
filterlarge.width=100;
%% Pixel
pixel = pixel3D('width',width);
largepixel = pixel3D('width',100);


%% Run simulation for each fnumber and chief ray angle

fig=figure(5);clf;  
fig.WindowState='maximized'
count=1;

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
        
        % Visualize
        subplot(numel(fnumbers),numel(cradegs),count)
        displaySetup3D(filter,'wavepacket',incident_wavepacket,'pixel',pixel);
        title(['f/' num2str(fnumbers(f)) ' - chief ray = ' num2str(cradegs(a)) ' deg'])
        pause(0.1);
        
        Twave(:,a,f)=transmittanceTiny3D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixel);
        lensvignet=lensVignetted(eo16.exitpupil,eo16.P,eo16.h);
        

        Tray(:,a,f) = transmittanceTinyRayFocused(lensvignet,nair,effective_index,nsub,fwhm2reflectance(normalized_fwhm),width,targetcwl,wavelengths,conedeg,cradeg,polarization,accuracy,true);


        

        wavepacketIdealLens =  wavepacket3DLens(conedeg,cradeg,azimuth_deg);
                Trayideal(:,a,f)=transmittanceTinyRayFocused(lensIdeal,nair,effective_index,nsub,fwhm2reflectance(normalized_fwhm),width,targetcwl,wavelengths,conedeg,cradeg,polarization,accuracy,true);
        Tinf(:,a,f)=transmittanceTinyRayFocused(lensvignet,nair,effective_index,nsub,fwhm2reflectance(normalized_fwhm),100,targetcwl,wavelengths,conedeg,cradeg,polarization,accuracy,true);
               
        count=count+1;
    end
end



%% Plot transmittance
cmap = hot;
s=size(cmap,1);
color{1}=cmap(1,:);
color{2}=cmap(round(0.4*s),:);
color{3}=cmap(round(0.6*s),:)
color{4}=cmap(round(0.65*s),:)

maxnorm = @(x) x/max(x)
count=1;

  fig=figure(f);clf;  hold on;
  fig.Position=[752 154 1021 466];  
for f =1:numel(fnumbers)
    
  

    for a=1:numel(cradegs)
        subplot(numel(fnumbers),numel(cradegs),count); hold on;
        hinf(a)=plot(wavelengths,maxnorm(Tinf(:,a,f)),':','color','k','linewidth',1)
        hwave(a)=plot(wavelengths,maxnorm(Twave(:,a,f)),'color','r','linewidth',1)
        hray(a)=plot(wavelengths,maxnorm(Tray(:,a,f)),'color','b','linewidth',1)
%        hrayideal(a)=plot(wavelengths,maxnorm(Trayideal(:,a,f)),'color','m','linewidth',1)
        
        
        ylabel('Transmittance')
        xlabel('Wavelength (µm)')
        xlim([0.72 0.82])
        ylim([0 1])
        
        title(['f/' num2str(fnumbers(f)) ' - chief ray = ' num2str(cradegs(a)) ' deg'])
        count=count+1;
    end
    

    legend([hwave(1) hray(1) hinf(1)],'Tiny Wave','Tiny Ray','Infinite filter','location','best')
end

