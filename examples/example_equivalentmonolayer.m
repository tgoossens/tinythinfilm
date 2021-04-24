%% Comparison of tiny Fabry-Pérot filter and equivalent tiny monolayer
% This example demonstrated how the transmittance of a tiny dielectric Fabry-Pérot filter
% differs from equivalent monolayer response.
%
% It shows that a well chosen equivalent monolayer approximates the response of a Fabry-pérot well. That is, given that the Fabry-Pérot is nicely approximated by a Lorentzian.
%
%
% Copyright Thomas Goossens

clear; close all;


addpath('/home/thomas/Documents/tinyfilters/research/wavepacket')

    
%% Choose simulation options

polarisation = 'unpolarized';
accuracy = 7; 
wavelengths=linspace(0.73,0.85,300); % µm
angles = [0 10 15 20 ]; 


%% Create dielectric Fabry Perot filter using two materials

% Target central wavelength (where the filter is centered)
targetcwl = 0.800; %micron


%% Refractive indeices of indicent medium and substrate
nair=1;
nsub=3.56; %silicon substarte


%% Tiny Fabry-Pérot filter design
nl = 1.4; % low refractive index
nh = 2.4 % high refractive index

dh = targetcwl/(4*nh);%quarterwave 
dl = targetcwl/(4*nl);%quarterwave 

n = [nh nl nh nl nh nl nh [nl nl] nh nl nh nl nh nl nh];
thickness = [dh dl dh dl dh dl dh [dl dl] dh dl dh dl dh dl dh];

width=5.5; %micron
filter=tinyfilter(nair,n,nsub,thickness,width);

%% Calculation of the ffective refractive index (cfr. MACLEOD)
neff=nl*sqrt(1/(1-nl/nh+nl^2/nh^2));


%% Equivalent filter
% Calculate normalized FWHM from infinite filter at normal incidence
Tclassic(:,1)=classictransmittance(filter,0,wavelengths,polarisation);
normalized_fwhm=fwhm(wavelengths,Tclassic(:,1))/targetcwl;
filter_equivalent=tinyfilter_equivalentmonolayer(targetcwl,normalized_fwhm,neff,width,nair,nsub);



%% Run simulation for each angle
for a=1:numel(angles)
    

    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    Tclassic(:,a)=classictransmittance(filter,angles(a),wavelengths,polarisation);


    Tmono(:,a)=tinytransmittance(filter_equivalent,angles(a),wavelengths,polarisation,accuracy);
    Ttiny(:,a)=tinytransmittance(filter,angles(a),wavelengths,polarisation,accuracy);

    
end


%% Plot transmittance

cmap = hot;
s=size(cmap,1);
color{1}=cmap(1,:);
color{2}=cmap(round(0.45*s),:);
color{3}=cmap(round(0.5*s),:);
color{4}=cmap(round(0.6*s),:);
color{5}=cmap(round(0.66*s),:)


figure(1);clf;  hold on;
for a=1:numel(angles)
    subplot(2,2,a); hold on;
    htiny=plot(wavelengths,Ttiny(:,a),'color',color{2},'linewidth',2)
    hmono=plot(wavelengths,Tmono(:,a),'-','color','k','linewidth',2)
    hclassic=plot(wavelengths,Tclassic(:,a),':','color','k','linewidth',1.5)
    ylabel('Transmittance')
    xlabel('Wavelength (µm)')
end




legh=legend([htiny hmono hclassic],'Tiny Fabry-Pérot','Equivelent monolayer','Infinite filter')
legh.Orientation='horizontal';
legh.Position=[0.1511 0.9469 0.7625 0.0429];

box on

























