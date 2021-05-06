%% Comparison of Tiny filter models with FDFD simulations for a thin-filmfilter sandwiched between two perfectly reflective boundaries
%
% The model including the reflective boundaries approximates the FDFD
% simulation well (especially the bimodality). The tiny wave optics model
% based on aperture diffraction does not apply to this case.


% Copyright Thomas Goossens

%%
clear; close all;



%% Load FDFD simulations
angles = [0 10 15 20];
FDFD={};
for a=1:numel(angles)
    FDFD =load(['./data/reflective_power_5500_' num2str(angles(a)) '.mat']); 
    Tfdfd(:,a)=FDFD.power.filter2.transmitted./FDFD.power.filter2.incident;
    wavelengths_fdfd=FDFD.wls /1000;%to micron
end

%% Create dielectric Fabry Perot filter using two materials

% Target central wavelength
targetcwl = 0.720; %micron

nair=1;
nsub=3.67; %silicon substarte

nl = 1.5; % low refractive index
nh = 2.4; % high refractive index

dh = targetcwl/(4*nh);%quarterwave
dl = targetcwl/(4*nl);%quarterwave

extra=0.1;
n = [nh nl nh nl nh nl nh nl nh (nl) nl   nh nl nh nl nh nl nh nl nh];
thickness = [dh dl dh dl dh dl dh dl dh (dl) dl   dh dl dh dl dh dl dh dl dh];

neff=nl/sqrt((1-nl/nh+nl^2./nh^2));
width=5.5 %micron



%% Define simulation
filter=tinyfilterCreate(nair,n,nsub,thickness,width);



%% Choose simulation options

polarization = 's';

accuracy = 8;
wavelengths=linspace(0.66,0.75,300); % µm


%% Run simulation for each angle
for a=1:numel(angles)
    
    
    
    
    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    Ttinyrefl(:,a)=transmittanceTiny2DReflectiveBoundaries(filter,angles(a),wavelengths,polarization);
    Ttiny(:,a)=transmittanceTiny2DCollimated(filter,angles(a),wavelengths,polarization,accuracy);
    Tinf(:,a)=transmittanceInfinite(filter,angles(a),wavelengths,polarization);
    
    
end


%% Plot transmittance
% There there is a drop in transmittance and an increase in FWHM
%



fig=figure(1);clf;  hold on;
fig.Position=[385 355 1215 383];
for a=1:numel(angles)
    subplot(2,2,a); hold on;
    hfdfd(a)=plot(wavelengths_fdfd,Tfdfd(:,a),'.-','color','k','linewidth',1,'markersize',10) 
    htiny(a)=plot(wavelengths,Ttiny(:,a),'color',[1 0.8 0.5],'linewidth',2)    
    htinyrefl(a)=plot(wavelengths,Ttinyrefl(:,a),'color','r','linewidth',2) 
    hinf(a)=plot(wavelengths,Tinf(:,a),':','color','k','linewidth',1.5)
    
    ylabel('Transmittance')
    xlabel('Wavelength (µm)')
    title([num2str(angles(a)) ' deg'])
    box on
end
legend([ hfdfd(1) htinyrefl(1) htiny(1) hinf(1)],'Numerical (FDFD)','Tiny Reflective','Tiny Diffractive','Infinite filter','location','best')































