%% FDFD validation for the case of a filter spanning multiple pixels
% (see geometry.png)
%
% A single filter spans 4 pixels and has different filters.
%
% From the FDFD simulations it follows that each pixel has a transmittance
% with a different angular dependency. The tiny filter model approximatites
% the transmittances very well.
%
%
%
% Copyright Thomas Goossens

%%
clear; close all;



%% Load FDFD simulations
angles = [ 20];
extras = [50 20 10 5];
widths = [3000];%nm


for px=3:6
    i=px-2;
    for e=1:numel(extras)

        FDFD =load(['./data/power_' num2str(widths) '_' num2str(angles) '_extra' num2str(extras(e)) '.mat']);
        Tfdfd(:,e,i)=FDFD.power.pixel{px}.transmitted./FDFD.power.pixel{px}.incident;
        Tfdfdneighbour(:,e)=FDFD.power.pixel{2}.transmitted./FDFD.power.pixel{2}.incident;
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


n = [nh nl nh nl nh nl nh nl nh (nl) nl   nh nl nh nl nh nl nh nl nh];
thickness = @(extra) [dh dl dh dl dh dl dh dl dh (dl-extra) dl   dh dl dh dl dh dl dh dl dh];

neff=nl/sqrt((1-nl/nh+nl^2./nh^2));






%% Choose simulation options

% FDFD simulation was for s-polarized light
polarization = 's';

accuracy = 7;
wavelengths=linspace(0.6,0.75,300); % µm


%% Options

filterwidths=[5*3 3];
pixelrange=0.5*[-filterwidths(1) filterwidths(1)]

%% Run simulation for each angle
% Loop over the four pixels and simulate the transmittance for 

fig=figure(5);clf;
fig.WindowState='maximized'
count=1;
for e=1:numel(extras)
    for px=3:6

        
        width=widths/1000; %nm->micron
        filterwidth=4*width;
        filter=tinyfilterCreate(nair,n,nsub,thickness(0),filterwidth);
%        filter=tinyfilterCreateEquivalent(0.720,0.007,1.71,5.5,1,nsub);
        filterneighbor=tinyfilterCreate(nair,n,nsub,thickness(extras(e)/1000),3*width);
            
        %% Determine equivalent model
        cwls_micron = [713.8 720]/1000;
        %cwls_micron = [708 720]/1000;
        R=fwhm2reflectance(0.007)
        
        Tray(:,e,px-2)=transmittanceTinyRayEquivalent_crosstalk(nair,neff,nsub,R,filterwidths,cwls_micron,wavelengths,angles,polarization,accuracy,pixelrange);
        
        
        %% Simulate
        disp(['Simulate tiny filter: width ' num2str(width) ' µm - ' num2str(extras(e)) ' extra']);
        
        
        
        
        
        % Define incident light
        wavepacket=wavepacket2DCollimated(angles,nair);
        
        % Define pixel kernel
        i=px-3;
        pixel=pixel2D('range',[-2*width+i*width, -width+i*width]);
        pixelneighbour=pixel2D('range',filterneighbor.width/2+[-width 0]);

        
        % Visualize
        subplot(numel(extras),4,count);
        displaySetup2D(filter,'wavepacket',wavepacket,'pixel',pixel);
        title(['Pixel ' num2str(px) ' - ' num2str(angles) ' deg' ])
        pause(0.1);
        
        % Simulate Tiny Filter
        Ttiny(:,e,px-2)=transmittanceTiny2D(filter,wavepacket,wavelengths,polarization,accuracy,pixel);
        Ttiny_neighbor(:,e,px-2)=transmittanceTiny2D(filterneighbor,wavepacket,wavelengths,polarization,accuracy,pixelneighbour);
        Tinf(:,e)=transmittanceInfinite(filter,angles,wavelengths,polarization);
        
        
        count=count+1;
    end
end


%% Plot transmittance
fig=figure(1);clf;  hold on;
count=1;
for e=1:numel(extras)
    for px=3:6
        i=px-2;
        fig.Position= [468 215 1105 731];
        
        subplot(numel(extras),4,count); hold on;
        htiny(e)=plot(wavelengths,Ttiny(:,e,i),'color',[1 0.1 0.5],'linewidth',2)
        hraycross(e)=plot(wavelengths,0.2*Tray(:,e,i),'color',[0.1 0.6 1],'linewidth',2)
%        htinye(e)=plot(wavelengths,Ttiny_neighbor(:,e,i),'-','color','m','linewidth',2)
        hfdfd(e)=plot(wavelengths_fdfd,Tfdfd(:,e,i),'.-','color','k','linewidth',1,'markersize',8)
        hinf(e)=plot(wavelengths,Tinf(:,e),':','color','k','linewidth',1.5)
        %plot(wavelengths_fdfd,Tfdfdneighbour(:,e),'g--')
        ylabel(num2str(extras(e)))
        xlabel('Wavelength (µm)')
        if(e==1);
            title(['Pixel ' num2str(px)])
        end
        box on
        
        count=count+1;
        ylim([0 1])
    end
    
    legh=legend([hfdfd(1) htiny(1) hinf(1)],'FDFD simulation','Tiny Wave model','Infinite filter')
    legh.Orientation='horizontal'
end





























