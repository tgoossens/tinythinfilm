%% Comparison of Tiny filter models with FDFD simulations for a thin-film filter sandwiched two other thin-film filters 
%
%  Only the the transmittance of the central filter (filter2) is analyzed
%


% Copyright Thomas Goossens

%%
clear; close all;



%% Load FDFD simulations
angles = [0 10 15 20];

widths = [5500 2000];%nm

for w=1:numel(widths)
    for a=1:numel(angles)
        FDFD =load(['./data/reflective_power_' num2str(widths(w)) '_' num2str(angles(a)) '.mat']);
        Tfdfd(:,a,w)=FDFD.power.filter2.transmitted./FDFD.power.filter2.incident;
        wavelengths_fdfd=FDFD.wls /1000;%to micron
    end
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






%% Choose simulation options

polarization = 's';

accuracy = 8;
wavelengths=linspace(0.66,0.75,300); % µm


%% Run simulation for each angle
for w=1:numel(widths)
    for a=1:numel(angles)
    
    width=widths(w)/1000; %nm->micron
    filter=tinyfilterCreate(nair,n,nsub,thickness,width);   
    
    %% Simulate
    disp(['Simulate tiny filter: width ' num2str(width) ' µm -' num2str(angles(a)) ' deg']);
    
    Ttinyrefl(:,a,w)=transmittanceTiny2DReflectiveBoundaries(filter,angles(a),wavelengths,polarization);
    Tinf(:,a,w)=transmittanceInfinite(filter,angles(a),wavelengths,polarization);
    
    
end
end


%% Plot transmittance
% There there is a drop in transmittance and an increase in FWHM
%




for w=1:numel(widths)
    fig=figure(w);clf;  hold on;
    fig.Position=[385 355 1215 383];
    for a=1:numel(angles)
        subplot(2,2,a); hold on;
        htiny(a)=plot(wavelengths,Ttinyrefl(:,a,w),'color',[1 0.2 0.2],'linewidth',2)
        hfdfd(a)=plot(wavelengths_fdfd,Tfdfd(:,a,w),'.-','color','k','linewidth',1,'markersize',10)
        
        hinf(a)=plot(wavelengths,Tinf(:,a,w),':','color','k','linewidth',1.5)
        
        ylabel('Transmittance')
        xlabel('Wavelength (µm)')
        title([num2str(angles(a)) ' deg'])
        box on
    end
    legend([ hfdfd(1) htiny(1) hinf(1)],'Numerical (FDFD)','Tiny Diffractive','Infinite filter','location','best')
    
end





























