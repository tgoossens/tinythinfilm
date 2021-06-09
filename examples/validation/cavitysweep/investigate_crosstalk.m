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
angles =  [ 5 10 20];
    
extras = [10];
widths = [3000];%nm


for a = 1:numel(angles)
for px=3
    i=px-2;
    for e=1:numel(extras)
        FDFD =load(['./data/power_' num2str(widths) '_' num2str(angles(a)) '_extra' num2str(extras(e)) '.mat']);
        Tfdfd{a}=FDFD.power.pixel{px}.transmitted./FDFD.power.pixel{px}.incident;
        Tfdfdneighbour{a}=FDFD.power.pixel{2}.transmitted./FDFD.power.pixel{2}.incident;
        wavelengths_fdfd{a}=FDFD.wls(1:numel(Tfdfdneighbour{a})) /1000;%to micron
    end
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

filterwidths=[3*2 3];
pixelrange=0.5*[-filterwidths(2) filterwidths(2)]



%% Run simulation for each angle
% Loop over the four pixels and simulate the transmittance for 

fig=figure(5);clf;
fig.WindowState='maximized'
count=1;

for a=1:numel(angles)
    angledeg=angles(a);
    
for e=1
    for px=3
        
        
        width=widths/1000; %nm->micron
        filterwidth=4*width;
        filter=tinyfilterCreate(nair,n,nsub,thickness(0),filterwidth);
        %        filter=tinyfilterCreateEquivalent(0.720,0.007,1.71,5.5,1,nsub);
        filterneighbor=tinyfilterCreate(nair,n,nsub,thickness(extras(e)/1000),3*width);
        
        %% Determine equivalent model
        

        Tn0=transmittanceInfinite(filterneighbor,0,wavelengths,polarization);
        Tn20=transmittanceInfinite(filterneighbor,20,wavelengths,polarization);
        
        cwl_neighbour = wavelengths(find(max(Tn0)==Tn0))
        normalized_fwhm_neighbour = fwhm(wavelengths,Tn0)/cwl_neighbour

        
        cwls_micron = [cwl_neighbour 0.720]
        R=[fwhm2reflectance(normalized_fwhm_neighbour) fwhm2reflectance(0.007)]
        Tray(:,a)=transmittanceTinyRayEquivalent_core(nair,neff,nsub,R(2),filterwidths(2),cwls_micron(2),wavelengths,angledeg,polarization,accuracy,pixelrange,true);
        Tcross(:,a)=transmittanceTinyRayEquivalent_crosstalk(nair,neff,nsub,R,filterwidths,cwls_micron,wavelengths,angledeg,polarization,accuracy,pixelrange,true);
        Tcrossincoherent(:,a)=transmittanceTinyRayEquivalent_crosstalk(nair,neff,nsub,R,filterwidths,cwls_micron,wavelengths,angledeg,polarization,accuracy,pixelrange,false);
        
        
        %% Simulate
        disp(['Simulate tiny filter: width ' num2str(width) ' µm - ' num2str(extras(e)) ' extra']);
        
        
        
        
        
        % Define incident light
        wavepacket=wavepacket2DCollimated(angledeg,nair);
        
        % Define pixel kernel
        i=px-3;
        pixel=pixel2D('range',[-2*width+i*width, -width+i*width]);
        %pixelneighbour=pixel2D('range',filterneighbor.width/2+[-width 0]);

        
        % Visualize
        subplot(numel(extras),4,count);
        displaySetup2D(filter,'wavepacket',wavepacket,'pixel',pixel);
        title(['Pixel ' num2str(px) ' - ' num2str(angledeg) ' deg' ])
        pause(0.1);
        
        % Simulate Tiny Filter
        Ttiny(:,a)=transmittanceTiny2D(filter,wavepacket,wavelengths,polarization,accuracy,pixel);
        %Ttiny_neighbor(:,a)=transmittanceTiny2D(filterneighbor,wavepacket,wavelengths,polarization,accuracy,pixelneighbour);
        Tinf(:,a)=transmittanceInfinite(filter,angledeg,wavelengths,polarization);
        Tinfneighbour(:,a)=transmittanceInfinite(filterneighbor,angledeg,wavelengths,polarization);
        
        
        count=count+1;
    end
end
end

%% Plot transmittance

% WHY IS THERE A 20% additional DROP?
unknownscalar=1



maxnorm = @(x)x*diag(max(x).^-1);
maxnorm =@(x)x;
fig=figure(1);clf;  hold on;
count=1;

for a=1:numel(angles)
for e=1
    for px=3
        i=px-2;
        fig.Position= [468 215 1105 731];
        fig.Position=[264 255 998 392];
        subplot(1,numel(angles),count); hold on;
        htiny(e)=plot(wavelengths,maxnorm(Ttiny(:,a)),':','color',[1 0.1 0.5],'linewidth',1)
        hraycross(e)=plot(wavelengths,maxnorm(Tcross(:,a)),'color',[0.1 0.6 1],'linewidth',2)
        hraycross(e)=plot(wavelengths,maxnorm(Tcrossincoherent(:,a)),'--','color',[0.1 0.6 1],'linewidth',2)
%        hray(e)=plot(wavelengths,maxnorm(Tray(:,a)),'color','r','linewidth',2)
        hfdfd(e)=plot(wavelengths_fdfd{a},maxnorm(Tfdfd{a}),'.-','color','k','linewidth',1,'markersize',8)
        hinf(e)=plot(wavelengths,Tinf(:,a),':','color','k','linewidth',1.5)
        hinfneigh(e)=plot(wavelengths,Tinfneighbour(:,a),':','color','k','linewidth',1.5)
        %plot(wavelengths_fdfd,Tfdfdneighbour(:,a),'g--')
        ylabel(num2str(extras(e)))
        xlabel('Wavelength (µm)')
        title([num2str(angles(a)) ' deg'])

        box on
        
        count=count+1;
        ylim([0 1])
    end
    

end

end
%    legh=legend([hfdfd(1) htiny(1) hinf(1)],'FDFD simulation','Tiny Wave model','Infinite filter')
    %legh.Orientation='horizontal'












return
%% Check equivalent model


Tn0=transmittanceInfinite(filterneighbor,0,wavelengths,polarization);
Tn20=transmittanceInfinite(filterneighbor,20,wavelengths,polarization);

cwl_neighbour = wavelengths(find(max(Tn0)==Tn0))
normalized_fwhm = fwhm(wavelengths,Tn0)/cwl_neighbour

equivfilter=tinyfilterCreateEquivalent(cwl_neighbour,1.6*normalized_fwhm,1.71,3,1,nsub);

Tequiv = transmittanceInfinite(equivfilter,0,wavelengths,polarization)
Tequiv20 = transmittanceInfinite(equivfilter,20,wavelengths,polarization)
figure(5);clf;hold on;
plot(wavelengths,Tn0,'k');
plot(wavelengths,Tequiv,'r');

plot(wavelengths,Tn20,'k');
plot(wavelengths,Tequiv20,'r');






