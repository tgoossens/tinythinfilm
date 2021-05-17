%% Example of a Tiny Fabry-Pérot filter with fully reflective boundaries
% This example demonstrates how the transmittance of a tiny dielectric
% Fabry-Pérot filter with perfectly reflective boundaries
% differs from the infinitely wide filter response.
%
% The result shows that at larger angles multiple peaks occur because
% multiple eigenmodes exist. This corresponds to a discretization of the
% spatial frequency domain .
%
% Copyright Thomas Goossens

%%

clear; close all;
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

width=5.5 %micron



%% Define simulation
filter=tinyfilterCreate(nair,n,nsub,thickness,width);



%% Choose simulation options

polarization = 'unpolarized';

accuracy = 8;
wavelengths=linspace(0.66,0.75,300); % µm
angles = [0 10 20 30];

%% Run simulation for each angle


for a=1:numel(angles)

    
    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);

    [Ttinyrefl(:,a)]=transmittanceTiny2DReflectiveBoundaries(filter,angles(a),wavelengths,polarization);
    [Ttiny(:,a)]=transmittanceTiny2DCollimated(filter,angles(a),wavelengths,polarization,accuracy);
    Tinf(:,a)=transmittanceInfinite(filter,angles(a),wavelengths,polarization);
    
    
end


%% Plot transmittance
fig=figure(1);clf;  hold on;
fig.Position=[385 355 1215 383];
for a=1:numel(angles)
    subplot(2,2,a); hold on;
    hclassic=plot(wavelengths,Tinf(:,a),':','color','k','linewidth',1.5)
    htinyrefl(a)=plot(wavelengths,Ttinyrefl(:,a),'color','r','linewidth',2)
    htiny(a)=plot(wavelengths,Ttiny(:,a),'color','k','linewidth',2)
    
    ylabel('Transmittance')
    xlabel('Wavelength (µm)')
    title([num2str(angles(a)) ' deg'])
    box on
end
legend([htiny(1) htinyrefl(1) hclassic(1)],'Tiny Diffraction','Tiny Reflective','Infinite filter')































