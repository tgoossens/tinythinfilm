%% Filterbank WAVE and RAY simulation for collimated (high fnumber) light.
% 
%
% This script simulates the effect illuminating a filter array with
% collimated light and the corresponding effect on the measurement of a
% reflectance spectrum with a spectral camera.
%
% The filter design here is based on the IMEC SWIR Mosaic sensor (16 wavebands)
% (May 2021)
%
% A comparison is made between the Tiny RAY and WAVE optics model.
% The ray model gives almost instant results. The wave optics model takes longerer to compute
% This script shows that for the given sensor dimensions
%
% A user can check whether for his particular system parameters the wave
% and ray optics match well. The ray optics can then be used for very fast
% computation as there is also an explicit analytical description.
%
% Copyright Thomas Goossens




clear;


%%%%%%%%%% 1. SIMULATION SETUP %%%%%%%%%



%% Define filters
nbFilters = 1;
nbWavelengths=2^9;

wavelengths=linspace(1000,1700,nbWavelengths);

cwl = linspace(1100,1650,nbFilters);
nFWHM=10/1100
width=15000; %nm


filterwidth=width;

neff=1.7;
nsub=1;
nair=1;



% Simulation parameters
polarization = 's'
accuracy=9

angledeg = [0 10 15 20];
fnumber=6;
coneangle_deg = atand(1./(2*fnumber));
alpha_deg=0; % azimuth

%%
%%%%%%%%%% 2. SIMULATION OF TRANSMITTANCES %%%%%%%%%

tic
for a=1:numel(angledeg)
    
    for c=1:nbFilters
        disp(['CRA ' num2str(angledeg(a)) ' deg - Filter ' num2str(c)]),

        % Define equivalent monolayer filter
        filter=tinyfilterCreateEquivalent(cwl(c),nFWHM,neff,width,nair,nsub);
        
        % Define incident light
        wavepacket2d = wavepacket2DCollimated(angledeg(a),nair,filterwidth);

        % Calculate transmittance for the tiny wave, tiny ray and infinite
        % filter case
        
        Trayanalytical(:,c,a)=transmittanceTinyRayEquivalent(nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,angledeg(a),polarization,accuracy,true);
        Tray(:,c,a)=transmittanceTinyRayEquivalent(nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,angledeg(a),polarization,accuracy,false);  
                
                
        Tinf(:,c,a)=transmittanceInfinite(filter,angledeg(a),wavelengths,polarization);
    end
    
end
toc


%% Compare WAVE and RAY optics prediction of  filter transmittances 
% The ray and wave optics models match each other very well.

figure(2);clf;hold on;
for a=1:numel(angledeg)
   subplot(numel(angledeg),1,a); hold on;
   plot(wavelengths,Tray(:,:,a),'k')
   plot(wavelengths,Trayanalytical(:,:,a),'r.')    
   
       title(['f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg'])
       ylim([0 1])

    xlabel('Wavelength (nm)')
    ylabel('Transmittance(%)')
end