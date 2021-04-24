%% Example of a Tiny Fabry-Pérot filter
% This example demonstrated how the transmittance of a tiny dielectric Fabry-Pérot filter
% differs from the infinitely wide filter response.
%
% Copyright Thomas Goossens

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

width=5.5; %micron



%% Define simualtion
filter=tinyfilter(nair,n,nsub,thickness,width);



%% Choose simulation options

polarisation = 'unpolarized';

accuracy = 7;
wavelengths=linspace(0.73,0.85,300); % µm
angles = [0 5 10 15 20 ]; 

%% Run simulation for each angle
for a=1:numel(angles)
    
   
    
    
    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    [Ttiny]=tinytransmittance(filter,angles(a),wavelengths,polarisation,accuracy);
    T(:,a)=Ttiny;

    Tclassic(:,a)=classictransmittance(filter,angles(a),wavelengths,polarisation);

end


%% Plot transmittance
% There there is a drop in transmittance and an increase in FWHM
% 

cmap = hot;
s=size(cmap,1);
color{1}=cmap(1,:);
color{2}=cmap(round(0.45*s),:);
color{3}=cmap(round(0.5*s),:);
color{4}=cmap(round(0.6*s),:);
color{5}=cmap(round(0.66*s),:)


figure(1);clf;  hold on;
for a=1:numel(angles)
    htiny(a)=plot(wavelengths,T(:,a),'color',color{a},'linewidth',2)
    hclassic=plot(wavelengths,Tclassic(:,a),':','color',color{a},'linewidth',1.5)
end




%% Labeling

text(0.73,0.6092,sprintf('Transmittance for\ninfinitely wide filter'))
line( [0.7617   0.7790],[ 0.5684    0.4888],'color','k')

text(0.73,0.2092,sprintf('Transmittance for\ntiny filter'))
line( [ 0.7608    0.7729],[  0.1969    0.1255],'color','k')

legend([htiny],'0^\circ','5^\circ','10^\circ','15^\circ','20^\circ')



ylabel('Transmittance')
xlabel('Wavelength (µm)')
title('Tiny vs. infinite transmittance')
box on

























