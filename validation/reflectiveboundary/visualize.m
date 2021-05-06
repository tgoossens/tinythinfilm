clear;close all;


%% Load FDFD simulations
width=2; %micron

addpath(genpath('/home/thomas/Documents/tinyfilters/git/tinythinfilm/core'))

angles = []
T={}
try     T{end+1} = load('./reflective_power_2000_0.mat');     angles(end+1)=0;catch 

end
try     T{end+1} = load('./reflective_power_2000_10.mat');     angles(end+1)=10;catch 
end
try     T{end+1} = load('./reflective_power_2000_15.mat');    angles(end+1)=15; catch 
end
try     T{end+1} = load('./reflective_power_2000_20.mat');     angles(end+1)=20;catch 
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



filter=tinyfilterCreate(nair,n,nsub,thickness,width);

%% Choose simulation options

polarisation = 's';

accuracy = 9;
wavelengths=[660:1:750]*1e-3;



%% Run simulation for each angle
for a=1:numel(angles)
    

    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    Ttiny(:,a)=transmittanceTiny2DCollimated(filter,angles(a),wavelengths,polarisation,accuracy);
    Ttinyrefl(:,a)=tinytransmittance_reflecting_core(filter,angles(a),wavelengths,polarisation,accuracy);


    Tclassic(:,a)=transmittanceInfinite(filter,angles(a),wavelengths,polarisation);

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


fig=figure;clf;  hold on;
fig.Position=[437 64 1046 882];
for a=1:numel(T)
    subplot(numel(angles),1,a); hold on;

    hclassic=plot(wavelengths,Tclassic(:,a),':','color','k','linewidth',1)
    htiny(a)=plot(wavelengths,Ttiny(:,a),'color',color{5},'linewidth',2)    
    htinyrefl(a)=plot(wavelengths,Ttinyrefl(:,a),'-','color',color{2},'linewidth',2)    
    % FDFD
    F = T{a};
    hfdfd=plot(1e-3*F.wls(1:F.i),F.power.filter2.transmitted(1:F.i)./F.power.filter2.incident(1:F.i),'k.-','markersize',10)


    title([num2str(angles(a)) ' deg'])
end




%% Labeling

legend([hfdfd(1) hclassic(1) htiny(1) htinyrefl(1)],'FDFD','Infinite','Tiny','Tiny Reflective Boundaries')


ylabel('Transmittance')
xlabel('Wavelength (µm)')

box on



























