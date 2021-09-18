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
close all

%%%%%%%%%% 1. SIMULATION SETUP %%%%%%%%%

%%
addpath(genpath(['/home/thomas/Documents/imec/libs_hsimatlab-15-06-' ...
                 '2018']))

sc=cSensorCalibration('/home/thomas/Documents/imec/phd/references/calibrations/CMV2K-SSM5x5-600_1000-5.3.12.8.xml');

M = sc.AlgorithmHyperspectral;

cwlorig = sc.GetCentralWavelengths
B = sc.GetBandResponses';
[M0,cwl0] = M.GetProperties;



% Define filters
nbFilters = 25;
nbWavelengths=2^7;

wavelengths=linspace(600,950,nbWavelengths);

%cwl = linspace(650,950,nbFilters);
[cwl,index]=sort(cwlorig);
B=B(:,index);


width=5500; %nm
widths=[linspace(6000,2000,50) 2000 2000 2000 2000 2000 2000 ];


angledeg=0;

neff=1.7;
nsub=1;
nair=1;



% Simulation parameters
polarization = 'unpolarized';
fastapproximation=true;
accuracy=9

angledeg = [0 15];

fnumber=2.8
coneangle_deg = atand(1./(2*fnumber));
alpha_deg=0; % azimuth

%% Reflectance of sample
% We use a sinusoidal reflectance.

clear R;
cyclespernm=0.01;
R = 1+sin(2*pi*wavelengths*cyclespernm)';
%S=csvread('samples/erbium_oxide.csv')
S=load('samples/FGB67.txt')

R=interp1(S(:,1),S(:,2),wavelengths');


%% Video prep

% create the video writer with 1 fps
for a=1:numel(angledeg)
    writerObj{a} = VideoWriter(['animation-' num2str(angledeg(a)) 'deg.avi']);
    writerObj{a}.FrameRate = 10;
    % open the video writer
    open(writerObj{a});
    
    figure(a); clf; hold on;
    set(gcf,'color','w');
    fig.Position=[376 153 658 352];
        
end
tic

%%

for wi=1:numel(widths)
    width=widths(wi);
     
    
    %%
    %%%%%%%%%% 2. SIMULATION OF TRANSMITTANCES %%%%%%%%%
    eo16=load('lens-eo16.mat').lens;
    
    %figure(1);clf;hold on;
    wlmeasured = 400:1000; wlsel=300:numel(wlmeasured);
    for a=1:numel(angledeg)
        for c=1:nbFilters
            CWL(c) = wlmeasured(find(B(:,c)==max(B(wlsel,c))));
            %nFWHM = fwhm(wlmeasured(wlsel),B(250:end,c))/CWL(c);
            nFWHM=0.3/100;
            %lens=lensVignetted(eo16.exitpupil,eo16.P,eo16.h);
            lens=lensIdeal();
            Trayvignet(:,c,a)=real(transmittanceTinyRayFocused(lens,nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation));
            largewidth=50*width;
            
            Tinf(:,c,a)=transmittanceTinyRayFocused(lens,nair,neff,nsub,1-pi*nFWHM,largewidth,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation);
        end
        
        
        
        
        for c=1:nbFilters
            tic
            disp(['CRA ' num2str(angledeg(a)) ' deg - Filter ' num2str(c)]),
            % Define equivalent monolayer filter
            filter=tinyfilterCreateEquivalent(cwl(c),nFWHM,neff,width,nair,nsub);
            
            % Calculate transmittance for the tiny wave, tiny ray and infinite
            % filter case
            
            
            
            % Wave optics
            accuracy_wave=7;
            wavepacket_lens = wavepacket3DLens(coneangle_deg,angledeg(a),0);
            %wavepacket_lens = wavepacket3DLensVignetted(coneangle_deg,angledeg(a),0,eo16.P,eo16.h,eo16.exitpupil);
            pixel=pixel3D('width',width);
            %            Twave(:,c,a)=transmittanceTiny3D(filter,wavepacket_lens,wavelengths,polarization,accuracy_wave,pixel);
            toc
        end
        
    end
    
    
    
    
    
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
        
        Brayvignet = Trayvignet(:,:,a);
        %   Bwave = Twave(:,:,a);
        Binf = Tinf(:,:,a);
        
        % Apply Flat Fielding
        drayvignet(:,a)= (Brayvignet'*R)  ./ (Brayvignet'*white);
        %dwave(:,a)= (Bwave'*R)  ./ (Bwave'*white);
        dinf(:,a)= (Binf'*R)./(Binf'*white);
    end
    
    
    
    
    %% Plot Reflectance values
    % It observed that, for the chosen spectrum, there is not much different
    % between the measured reflectance predicted by the tiny filter vs infinite
    % filter model.
    
    
    
    
    
    for a=1:numel(angledeg)
       
        figure(a); hold on;
        cra=deg2rad(angledeg(a));
        cone=deg2rad(coneangle_deg);
        cwlcorrect= cwl*(1-cra.^2/(2*neff^2) -cone.^2/(4*neff^2));
        htiny=plot(cwlcorrect,drayvignet(:,a),'.-','color',[1 0.2 0.2],'markersize',14,'linewidth',2)
        %hwave(a)= plot(cwlcorrect,dwave(:,a),'b.-','markersize',14,'linewidth',1.5)
        hinf= plot(cwlcorrect,dinf(:,a),'.-','markersize',14,'linewidth',1,'color',0.2 *[1 1 1])
        %href(a)=plot(wavelengths,R,'k:')
        
        % Labels
        %title(['\textbf{f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg}'])
        xlabel('Wavelength (nm)')
        ylabel('Transmittance(%)')
    
        
        
        
        xlim([650 900])
        
        
      
        
        %subplot(1,2,1)
        
        
        legh=legend([ htiny(1) hinf(1)],'Tiny Filter','Infinite Filter')
        %legh.Orientation='horizontal';
        legh.Location='north';
        legh.Box='off'
        legh.FontSize=12
        legh.Interpreter='latex'
        legh.Position=[0.4267 0.6717 0.2211 0.1219];
        
        
        
            
        % Pixel size label
        text(740,95,['\textbf{' sprintf('Pixel width = %.2f',widths(wi)/1000) ' $\mu$m}'])
        %text(740,95,[sprintf('Pixel width = %.2f',widths(wi)/1000)])
        
        
        
        
        set(findall(gcf,'-property','FontSize'),'FontSize',13)
        set(findall(gcf,'-property','interpreter'),'interpreter','latex')
        
        % Write video
     
        fig=gcf;
        fig.Position=[376 153 658 352];
        
        ax = gca;
   
        
        frame =   getframe(gcf);
        
        writeVideo(writerObj{a}, frame);
        clf;
        
    end
end
%%
for a=1:numel(angledeg)
    close(writerObj{a});
end
totaltime=toc





