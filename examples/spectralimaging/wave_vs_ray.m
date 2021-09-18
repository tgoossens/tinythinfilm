%% Wave vs ray model

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
widths=linspace(6000,2000,50);
widths=5500



neff=1.7;
nsub=3.67;
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

S=load('samples/FGB67.txt')

R=interp1(S(:,1),S(:,2),wavelengths');



%%

csel=1:nbFilters;

for wi=1:numel(widths)
    width=widths(wi);
    
    %%
    %%%%%%%%%% 2. SIMULATION OF TRANSMITTANCES %%%%%%%%%
    %eo16=load('lens-eo16.mat').lens;
    
    %figure(1);clf;hold on;
    wlmeasured = 400:1000; wlsel=300:numel(wlmeasured);
    for a=1:numel(angledeg)
        for c=csel
            CWL(c) = wlmeasured(find(B(:,c)==max(B(wlsel,c))));
            %nFWHM = fwhm(wlmeasured(wlsel),B(250:end,c))/CWL(c);
            nFWHM = 0.3/100;
            
            %lens=lensVignetted(eo16.exitpupil,eo16.P,eo16.h);
            lens=lensIdeal();
            Trayvignet(:,c,a)=real(transmittanceTinyRayFocused(lens,nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation));
            largewidth=50*width;
            
            Tinf(:,c,a)=transmittanceTinyRayFocused(lens,nair,neff,nsub,1-pi*nFWHM,largewidth,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation);
        end
        
        
        
        
        for c=csel
            tic
            disp(['CRA ' num2str(angledeg(a)) ' deg - Filter ' num2str(c)]),
            % Define equivalent monolayer filter
            filter=tinyfilterCreateEquivalent(cwl(c),nFWHM,neff,width,nair,nsub);
            
            % Calculate transmittance for the tiny wave, tiny ray and infinite
            % filter case
            
            
            
            % Wave optics
            accuracy_wave=8;
            wavepacket_lens = wavepacket3DLens(coneangle_deg,angledeg(a),0);
            %wavepacket_lens = wavepacket3DLensVignetted(coneangle_deg,angledeg(a),0,eo16.P,eo16.h,eo16.exitpupil);
            pixel=pixel3D('width',width);
            Twave(:,c,a)=transmittanceTiny3D(filter,wavepacket_lens,wavelengths,polarization,accuracy_wave,pixel);
            toc
        end
        
    end
    
end

%%
Trayvignet(isnan(Trayvignet))=0;
maxnorm = @(x) x*diag(max(x).^-1);
figure;
subplot(121)
plot(maxnorm(Twave(:,:,1)),'k');hold on;plot(maxnorm(Trayvignet(:,:,1)),'b');hold on;
subplot(122)
plot(maxnorm(Twave(:,:,2)),'k');hold on;plot(maxnorm(Trayvignet(:,:,2)),'b');hold on;
