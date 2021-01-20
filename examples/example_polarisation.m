% EXAMPLE POLARISATION
%
% This examples demonstates the different transmittance for s,p, and unpolarized light
% It is seen that the transmittance for  p polarized light collapses less than quickly than the s-polarized light.
% 
% Copyright (c) Thomas Goossens

clear; close all;



%% Create dielectric Fabry Perot filter using two materials

% Target central wavelength
targetcwl = 0.8; %micron


nair=1;
nsub=3.56; %silicon substarte

nl = 1.45; % low refractive index
nh = 2.4; % high refractive index

dh = targetcwl/(4*nh);%quarterwave 
dl = targetcwl/(4*nl);%quarterwave 

nstack = [nh nl nh nl nh nl nh  [nl nl]  nh nl nh nl nh nl nh];
thickness = [dh dl dh dl dh dl dh  [dl dl]  dh dl dh dl dh dl dh];

width=5.5 %micron


filter=tinyfilter(nair,nstack,nsub,thickness,width);



%% Choose simulation options


accuracy = 7;
wavelengths=linspace(0.7,0.85,300); % µm
angles = [0 10 20 30 40]; 


%% Run simulation for each angle
for a=1:numel(angles)
    

    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    [Ttinys]=tinytransmittance(filter,angles(a),wavelengths,'s',accuracy);
    [Ttinyp]=tinytransmittance(filter,angles(a),wavelengths,'p',accuracy);
    Ts(:,a)=Ttinys;    Tp(:,a)=Ttinyp;

    Tinf_s(:,a)=classictransmittance(filter,angles(a),wavelengths,'s');
    Tinf_p(:,a)=classictransmittance(filter,angles(a),wavelengths,'p');


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


fig=figure(1);clf; 
fig.Position=[677 428 800 308];
for a=1:numel(angles)
    subplot(131);hold on;
    title('S polarized')
    htiny(a)=plot(wavelengths,Ts(:,a),'color',color{a},'linewidth',2);
    hclassic=plot(wavelengths,Tinf_s(:,a),':','color',color{a},'linewidth',1.5);

    subplot(132);hold on;
    title('P polarized')
    plot(wavelengths,Tp(:,a),'color',color{a},'linewidth',2);
    plot(wavelengths,Tinf_p(:,a),':','color',color{a},'linewidth',1.5);

    
    subplot(133);hold on;
    title('Unpolarized (S+P)/2')
    plot(wavelengths,0.5*(Tp(:,a)+Ts(:,a)),'color',color{a},'linewidth',2);
    plot(wavelengths,0.5*(Tinf_p(:,a)+Tinf_s(:,a)),':','color',color{a},'linewidth',1.5);
    
end




%% Labeling
ylimits=[0 0.8];
xlimits=[0.7 0.82]
subplot(131)
legend([htiny],'0^\circ','10^\circ','20^\circ','30^\circ','40^\circ','location','best')
box on
ylabel('Transmittance')
ylim(ylimits)
xlim(xlimits)
subplot(132)
xlabel('Wavelength (µm)')
ylim(ylimits)
xlim(xlimits)
box on



subplot(133)
xlabel('Wavelength (µm)')
ylim(ylimits)
xlim(xlimits)
box on
























