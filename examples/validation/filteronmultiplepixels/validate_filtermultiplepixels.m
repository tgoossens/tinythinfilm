%% Comparison of Tiny filter models with FDFD simulations for a thin-film filter sandwiched two other thin-film filters 
%
%  Only the the transmittance of the central filter (filter2) is analyzed
%


% Copyright Thomas Goossens

%%
clear; close all;



%% Load FDFD simulations
angles = [0 10 15 20];

widths = [3000];%nm

w=1;

for px=3:6
    i=px-2;
    for a=1:numel(angles)

        FDFD =load(['./data/power_' num2str(widths(w)) '_' num2str(angles(a)) '.mat']);
        Tfdfd(:,a,i)=FDFD.power.pixel{px}.transmitted./FDFD.power.pixel{px}.incident;
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

% FDFD simulation was for s-polarized light
polarization = 's';

accuracy = 8;
wavelengths=linspace(0.66,0.75,300); % µm


%% Run simulation for each angle
for px=3:6
    for a=1:numel(angles)
    
        
    width=widths(w)/1000; %nm->micron
    filterwidth=4*width;
    filter=tinyfilterCreate(nair,n,nsub,thickness,filterwidth);   
    
    %% Simulate
    disp(['Simulate tiny filter: width ' num2str(width) ' µm - ' num2str(angles(a)) ' deg']);
    

    % Define incident light
    wavepacket=wavepacket2DCollimated(angles(a),nair);
    
    % Define pixel kernel
    i=px-3;
    pixel=pixel2D('range',[-2*width+i*width, -width+i*width]);

    
    % Simulate Tiny Filter
    Ttiny(:,a,px-2)=transmittanceTiny2D(filter,wavepacket,wavelengths,polarization,accuracy,pixel);
    Tinf(:,a)=transmittanceInfinite(filter,angles(a),wavelengths,polarization);
    
end
end


%% Plot transmittance
fig=figure(w);clf;  hold on;
count=1;
for a=1:numel(angles)
    for px=3:6
        i=px-2;
        fig.Position= [468 215 1105 731];
        
        subplot(4,4,count); hold on;
        htiny(a)=plot(wavelengths,Ttiny(:,a,i),'color',[1 0.1 0.5],'linewidth',2)
        hfdfd(a)=plot(wavelengths_fdfd,Tfdfd(:,a,i),'.-','color','k','linewidth',1,'markersize',8)
        hinf(a)=plot(wavelengths,Tinf(:,a),':','color','k','linewidth',1.5)
        
        ylabel('Transmittance')
        xlabel('Wavelength (µm)')
        if(a==1);
            title(['Pixel ' num2str(px)])
        end
        box on
        
        count=count+1;
        ylim([0 1])
    end
    
    legh=legend([hfdfd(1) htiny(1) hinf(1)],'FDFD simulation','Tiny Wave model','Infinite filter')
    legh.Orientation='horizontal'
end





























