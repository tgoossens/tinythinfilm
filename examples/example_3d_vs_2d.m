%% Example of difference between 2D and 3D simulation for collimated light at oblique incidence
%
% Copyright Thomas Goossens

clear; close all;


addpath('/home/thomas/Documents/tinyfilters/research/wavepacket')
addpath(genpath('../core'))

    
%% Choose simulation options

% Accuracy options
accuracy = 8;
hiaccuracy = 8; % for the large filter

% Light properties
wavelengths=linspace(0.73,0.85,100); % µm
polarization = 'unpolarized';

% Lens properties
angles=[0 5 10 20]; % chief ray angles

%% Define equivlent monolayer filter


% Target central wavelength
targetcwl = 0.800; %micron

% Incident medium and substrate
nair=1;
nsub=3.56; %silicon substarte

% Normalized bandwith
normalized_fwhm=1/100; % 1 percent

% Effective refractive index
effective_index=1.7;

% Filter width
filterwidth=5.5; %micron

% Equivalent filter create
filter=tinyfilterCreateEquivalent(targetcwl,normalized_fwhm,effective_index,filterwidth,nair,nsub);

%% Pixel
pixelkernel = pixel_fullwidth(filterwidth);


%% Run simulation for each fnumber and chief ray angle
    for a=1:numel(angles)
        angle=angles(a);
        

        %% Simulate
        
        
        
        disp(['Simulate tiny filter 2D collimated: angle = ' num2str(angle) ' deg']);
        
        Tinf(:,a) = transmittanceInfinite(filter,angle,wavelengths,polarization);
        
        wavepacket2D=  wavepacket2DCollimated(angle,nair,filterwidth);
        Ttiny2D(:,a)=transmittanceTiny2D(filter,wavepacket2D,wavelengths,polarization,accuracy,pixelkernel);
       
        disp(['Simulate tiny filter 3D collimated: angle = ' num2str(angle) ' deg']);
        azimuth_deg=0;
        wavepacket3D=  wavepacket3DCollimated(angle,azimuth_deg,nair,filterwidth);
        Ttiny3D(:,a)=transmittanceTiny3D(filter,wavepacket3D,wavelengths,polarization,accuracy,pixelkernel);
        


end



%% Plot transmittance
cmap = hot;
s=size(cmap,1);
color{1}=cmap(1,:);
color{2}=cmap(round(0.4*s),:);
color{3}=cmap(round(0.6*s),:)
color{4}=cmap(round(0.65*s),:)

maxnorm = @(x) x/max(x);
fig=figure(1);clf;  hold on;
fig.Position= [533 488 666 209];

    for a=1:numel(angles)
        subplot(1,numel(angles),a); hold on;
        hclassic(a)=plot(wavelengths,Tinf(:,a),':','color','k','linewidth',1)
        htiny2d(a)=plot(wavelengths,Ttiny2D(:,a),'color','k','linewidth',1)
        htiny3d(a)=plot(wavelengths,Ttiny3D(:,a),'color',color{2},'linewidth',1)
        
        
        ylabel('Transmittance')
        xlabel('Wavelength (µm)')
        xlim([0.72 0.82])
        ylim([0 1])
        title([num2str(angles(a)) ' deg'])
    end
    
    

    subplot
    legend([htiny2d(1) htiny3d(1)  hclassic(1)],'Tiny 2D','Tiny 3D', 'Infinite')


