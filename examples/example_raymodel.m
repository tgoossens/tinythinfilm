%% Example of a Tiny Fabry-Pérot filter calculated using a wave model and an equivalent ray optics model
%
% One can see that for the particular Fabry-Pérot design in this script,
% the ray model transmittance is almost indistinguishable from the
% wave-optics model.
%
% At normal incidence the ray optics model matches the infinite filter
% prediction perfectly. This is because diffraction is not included and
% hence, at normal incidence, an infinite number of reflection can occur.
% 
% Use this script to explore in what cases deviations become strong.
% Eg : filterwidth=2 %micron
%

%%

clear; close all;
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



neff=nl/sqrt((1-nl/nh+nl^2./nh^2));


% Filter width
filterwidth=5.5; %micron






%% Define simualtion
filter=tinyfilterCreate(nair,n,nsub,thickness,filterwidth);



%% Choose simulation options

polarisation = 'unpolarized';

accuracy = 6;
wavelengths=linspace(0.73,0.85,300); % µm
angles = [0 10 15 20 ]; 

%% Run simulation for each angle
for a=1:numel(angles)
    
   
    
    
    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    Ttiny(:,a)=transmittanceTiny2DCollimated(filter,angles(a),wavelengths,polarisation,accuracy);
    

    Tinf(:,a)=transmittanceInfinite(filter,angles(a),wavelengths,polarisation);
    
    % -- Ray model --
    %  1. Calculate equivalent monolayer parameters
    normalized_fwhm=fwhm(wavelengths,Tinf(:,1))/targetcwl;  %Normalized FWHM
    R=1-pi*normalized_fwhm; % Corresponding mirror reflectance
    % 2. Calculate the transmittance
    Tray(:,a)=transmittanceTinyRayEquivalent(nair,neff,nsub,R,filterwidth,targetcwl,wavelengths,angles(a)+eps,'s',accuracy,true);
        
end


%% Plot transmittance
fig=figure(1);clf;  hold on;
fig.Position=[385 355 1215 383];
for a=1:numel(angles)
    subplot(2,2,a); hold on;
    hinf=plot(wavelengths,Tinf(:,a),':','color','k','linewidth',1.5)
    htiny(a)=plot(wavelengths,Ttiny(:,a),'color','k','linewidth',2)
    hray(a)=plot(wavelengths,Tray(:,a),'-','color',[1 0.2 0.2],'linewidth',2)
    
    ylabel('Transmittance')
    xlabel('Wavelength (µm)')
    title([num2str(angles(a)) ' deg'])
    box on
end
legend([htiny(1) hray(1) hinf(1)],'Tiny Wave model','Tiny Ray Model','Infinite filter')
































