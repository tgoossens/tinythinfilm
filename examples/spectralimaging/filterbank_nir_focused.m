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
nbWavelengths=2^9;

wavelengths=linspace(600,950,nbWavelengths);

%cwl = linspace(650,950,nbFilters);
[cwl,index]=sort(cwlorig);
B=B(:,index);


%B=B(:,index);
nFWHM = 4.01/662
width=5500; %nm
width=2000
%width=5000

filterwidth=width;angledeg=0;

neff=1.7;
nsub=1;
nair=1;



% Simulation parameters
polarization = 'unpolarized';
fastapproximation=true;
accuracy=9

angledeg = [0 10];

fnumber=8
coneangle_deg = atand(1./(2*fnumber));
alpha_deg=0; % azimuth

%%
%%%%%%%%%% 2. SIMULATION OF TRANSMITTANCES %%%%%%%%%
eo16=load('lens-eo16.mat').lens;
tic
%figure(1);clf;hold on;
wlmeasured = 400:1000; wlsel=300:numel(wlmeasured);
for a=1:numel(angledeg)
        for c=1:nbFilters
            subplot(5,5,c); hold on;
            CWL(c) = wlmeasured(find(B(:,c)==max(B(wlsel,c))));
            nFWHM = fwhm(wlmeasured(wlsel),B(250:end,c))/CWL(c)
        Tray(:,c,a)=transmittanceTinyRayFocused(lensIdeal(),nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation);

        Trayvignet(:,c,a)=transmittanceTinyRayFocused(lensVignetted(eo16.exitpupil,eo16.P,eo16.h),nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation);
        
        largewidth=50*width;
eolens=lensVignetted(eo16.exitpupil,eo16.P,eo16.h);
        Tinf(:,c,a)=transmittanceTinyRayFocused(eolens,nair,neff,nsub,1-pi*nFWHM,largewidth,cwl(c),wavelengths,coneangle_deg,angledeg(a),polarization,accuracy,fastapproximation);
        
        
        %plot(400:1000,B(:,c)/max(B(:,c)),'k')
        %plot(wavelengths,Tinf(:,c,a),'r')
        %pause(0.1)
        end
        
      
    for c=1:nbFilters
        tic
        disp(['CRA ' num2str(angledeg(a)) ' deg - Filter ' num2str(c)]),
        % Define equivalent monolayer filter
        filter=tinyfilterCreateEquivalent(cwl(c),nFWHM,neff,width,nair,nsub);
        
        % Calculate transmittance for the tiny wave, tiny ray and infinite
        % filter case

        
          
        % Wave optics
        accuracy_wave=8;
%        wavepacket_lens = wavepacket3DLens(coneangle_deg,angledeg(a),0);
        wavepacket_lens = wavepacket3DLensVignetted(coneangle_deg,angledeg(a),0,eo16.P,eo16.h,eo16.exitpupil);
        pixel=pixel3D('width',filterwidth);
        Twave(:,c,a)=transmittanceTiny3D(filter,wavepacket_lens,wavelengths,polarization,accuracy_wave,pixel);
        toc
    end
    
end
toc


%% Compare WAVE and RAY optics prediction of  filter transmittances 
% The ray and wave optics models match each other very well.
maxnorm = @(x)x *diag(max(x).^-1);
Trayvignet(isnan(Trayvignet))=0;
maxnorm = @(x)x;
figure(2);clf;hold on;
for a=1:numel(angledeg)
   subplot(numel(angledeg),1,a); hold on;
   
%(wavelengths,maxnorm(Tray(:,:,a)),'k')
   plot(wavelengths,maxnorm(Trayvignet(:,:,a)),'color', [1 0.2 0.2])
   plot(wavelengths,maxnorm(Tinf(:,:,a)),'k:')
   
   plot(wavelengths,maxnorm(Twave(:,:,a)),'b')    
       title(['f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg'])
       ylim([0 1])

    xlabel('Wavelength (nm)')
    ylabel('Transmittance(%)')
end
%legend([hinf(1) hray(1) hwave(1)],'Infinite filter','Tiny Ray','Tiny Wave','location','north')

%%%%% REFLECTANCE MEASUREMENT %%%%%


%% Reflectance of sample
% We use a sinusoidal reflectance.

clear R;
cyclespernm=0.01;
R = 1+sin(2*pi*wavelengths*cyclespernm)';
%S=csvread('samples/erbium_oxide.csv')
S=load('samples/FGB67.txt')

R=interp1(S(:,1),S(:,2),wavelengths');
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
    Brayvignet = Trayvignet(:,:,a); 
%    Bwave = Twave(:,:,a);    
    Binf = Tinf(:,:,a); 
    
    % Apply Flat Fielding
    dray(:,a)= (Bray'*R)  ./ (Bray'*white);
    drayvignet(:,a)= (Brayvignet'*R)  ./ (Brayvignet'*white); 
    %dwave(:,a)= (Bwave'*R)  ./ (Bwave'*white);
    dinf(:,a)= (Binf'*R)./(Binf'*white);
end


%% Plot Reflectance values
% It observed that, for the chosen spectrum, there is not much different
% between the measured reflectance predicted by the tiny filter vs infinite
% filter model.


fig=figure(1);clf;hold on;
fig.Position=[610 87 988 721];
fig.Position=[610 87 510 535];
for a=1:2
    subplot(2,1,a); hold on;
    %plot(cwl,dray(:,a),'b.-')
    htiny(a)=plot(cwl,drayvignet(:,a),'.-','color',[1 0.2 0.2],'markersize',10)   
%   plot(cwl,dwave(:,a),'b.-')
%    plot(cwl,dinf(:,1),'k.-')
   hinf(a)= plot(cwl,dinf(:,a),'k.-','markersize',10)  
    href(a)=plot(wavelengths,R,'k:')
    
    title(['\textbf{f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg}'])
    xlabel('Wavelength (nm)')
    ylabel('Transmittance(%)')
    xlim([650 900])
end


set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

subplot(2,1,1)
legh=legend([href(1) htiny(1) hinf(1)],'Ground truth','Tiny Filter','Infinite Filter')
legh.Orientation='horizontal';
legh.Location='north';
legh.Box='off'
legh.FontSize=12
legh.Interpreter='latex'
legh.Position=[0.1282 0.9624 0.8234 0.0434];


saveas(gcf,['./fig/f' num2str(fnumber) '-width' num2str(width) '.png'])
%%
%saveas(gcf,['./fig/swir_collimated_width' num2str(width/1000) '-nfwhm' num2str(nFWHM*100) 'percent-brownglass.png'])



%%%%%%%%%% MEASUREMENT CHECK



%% Get files

fnum = [1.4 8];
for f = 1:numel(fnum)
    
    path=['/run/media/thomas/Passport/data/2019-02-18-mosaic-' ...
          'brownglass/eo16/f' num2str(fnum(f)) '/'];

    files = dir([path '/white/' '*.bmp']);
    white=0;

    for k = 1:length(files)
        white=((k-1)*white+single(imread(fullfile(files(k).folder,files(k).name))))/k;
    end

    % brown
    files = dir([path '/brown/' '*.bmp']);
    brown=0;

    for k = 1:length(files)
        brown=((k-1)*brown+single(imread(fullfile(files(k).folder,files(k).name))))/k;
    end

    T(:,:,:,f)=demosaic5x5(brown./white);

end






%%
addpath(genpath(['/home/thomas/Documents/imec/libs_hsimatlab_git/' ...
                 'libs_hsimatlab/lib_spectral_correction/angularity']))



exitpupil=21
for f = 1:numel(fnum)
 
Tcorr(:,:,:,f) = tmprod(T(:,:,:,f),M0,3);

custom_exit=1000*exitpupil; % large to avoid correction 

[Ttemp,cwlnew]= mosaic_correctcube(Tcorr(:,:,:,f),cwl0,0.001*atand(1/(2*fnum(f))),custom_exit,1.7,5.5e-3, ...
                    5,4,1);

Tnew(:,:,:,f)= Ttemp;
end

%%

%% COmpare measurements

f=2

fig=figure(10);clf;

fig.Position=[610 87 988 721];
mid = round(size(Tnew,1:2)/2);
cra = @(rowcol) atand(norm(rowcol-mid)*5*5.5e-3 / exitpupil);
% Center
rowcols=[100 200; 50 80; 20 20 ];

r=@(x) x:(x+20);

ref = mean(tens2mat(Tnew(r(mid(1)),r(mid(2)),:,f),3),2);
     
for i=1:size(rowcols,1)
    subplot(size(rowcols,1),1,i)
rowcol = rowcols(i,:);
plot(cwlnew,mean(tens2mat(Tnew(r(rowcol(1)),r(rowcol(2)),:,f),3),2),'r', ...
     'linewidth',2); hold on;
 plot(cwlnew,ref,'k');
 title(['f/' num2str(fnum(f)) '- ' num2str(cra(rowcol))])
 ylim([0 1])
 xlim([600 950])
end 





