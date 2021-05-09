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
nbFilters = 16;
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
polarization = 'unpolarized';
fastapproximation=true;
accuracy=9

angledeg = [0 10 15 20];
fnumber=2.8;
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
        Tray(:,c,a)=transmittanceTinyRayFocused(nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation);
        
        wavepacket_lens = wavepacket3DLens(coneangle_deg,angledeg(a),0,filter.width);
        %Twave(:,c,a)=transmittanceTiny3D(filter,wavepacket_lens,wavelengths,polarization,6,pixel_fullwidth(filter.width));
        Tinf(:,c,a)=transmittanceInfinite(filter,angledeg(a),wavelengths,polarization);
        toc
    end
    
end
toc


%% Compare WAVE and RAY optics prediction of  filter transmittances 
% The ray and wave optics models match each other very well.

figure(2);clf;hold on;
for a=1:numel(angledeg)
   subplot(numel(angledeg),1,a); hold on;
   plot(wavelengths,Tray(:,:,a),'k')
   plot(wavelengths,Twave(:,:,a),'r')    
       title(['f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg'])
       ylim([0 1])

    xlabel('Wavelength (nm)')
    ylabel('Transmittance(%)')
end


%%%%% REFLECTANCE MEASUREMENT %%%%%


%% Reflectance of sample
% We use a sinusoidal reflectance.

clear R;
cyclespernm=0.01;
R = 1+sin(2*pi*wavelengths*cyclespernm)';

%% Simulate flat fielding
% Flatfielding is a common procedure in spectral imaging to obtain
% reflectance values. It is obtained by measuring two images: an image of a
% 100% reflective (white) surface and the actual scene.
% The ratio between the two corrects for the unknown non-uniformities of
% the pixels and illumination.

% White reference 
white = ones(size(R));


% Calcualte reflectance values
for a=1:numel(angledeg)

    % Matrix containing filter responses for a specific angle
    Bray = Tray(:,:,a);
%    Bwave = Twave(:,:,a);    
    Binf = Tinf(:,:,a); 
    
    % Apply Flat Fielding
    dray(:,a)= (Bray'*R)  ./ (Bray'*white);
    %dwave(:,a)= (Bwave'*R)  ./ (Bwave'*white);
    dinf(:,a)= (Binf'*R)./(Binf'*white);
end



%% Plot Reflectance values
% It observed that, for the chosen spectrum, there is not much different
% between the measured reflectance predicted by the tiny filter vs infinite
% filter model.


fig=figure(1);clf;hold on;
fig.Position=[610 87 988 721];
for a=1:numel(angledeg)
    subplot(numel(angledeg),1,a); hold on;
    plot(cwl,dray(:,a),'r.-')
%    plot(cwl,dwave(:,a),'b.-')
    plot(cwl,dinf(:,a),'k.-')
    plot(wavelengths,R,'k:')
    
    title(['f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg'])
    xlabel('Wavelength (nm)')
    ylabel('Transmittance(%)')
    
end

legh=legend('Tiny Ray model','Tiny Wave model','Infinite filter','Ground truth')
legh.Orientation='horizontal';
legh.Location='north';
%%
%saveas(gcf,['./fig/swir_collimated_width' num2str(width/1000) '-nfwhm' num2str(nFWHM*100) 'percent-brownglass.png'])







