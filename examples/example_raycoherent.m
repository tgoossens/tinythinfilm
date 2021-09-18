%% Example of a Tiny Fabry-Pérot filter calculated using a wave model and an equivalent ray optics model
%
% One can see that for the particular Fabry-Pérot design in this script,
% the ray model transmittance is almost indistinguishable from the
% wave-optics model.
%
% At normal incidence the ray optics model matches the infinite filter
% prediction perfectly. This is because diffraction is not included and
% hence, at normal incidence, an infinite number of reflection can occur.
% 
% Use this script to explore in what cases deviations become strong.
% Eg : filterwidth=2 %micron
%

%%

clear; close all;
%% Create dielectric Fabry Perot filter using two materials

% Target central wavelength
targetcwl = 0.800; %micron
%targetcwl = 0.750; %micron

nair=1;
nsub=3.56; %silicon substarte

nl = 1.4; % low refractive index
nh = 2.9; % high refractive index

dh = targetcwl/(4*nh);%quarterwave 
dl = targetcwl/(4*nl);%quarterwave 

n = [nh nl nh nl nh nl nh [nl nl] nh nl nh nl nh nl nh];
thickness = [dh dl dh dl dh dl dh [dl dl] dh dl dh dl dh dl dh];



neff=nl/sqrt((1-nl/nh+nl^2./nh^2));


% Filter width
filterwidth=5.5; %micron
filterwidth=20





%% Define simualtion
filter=tinyfilterCreate(nair,n,nsub,thickness,filterwidth);



%% Choose simulation options

polarization = 'unpolarized';

accuracy = 8;
wavelengths=linspace(0.65,0.85,300); % µm
angles = [0 10 15 20 25 30]; 
%angles = [0 1 2 5 7 10]; 


%% Run simulation for each angle
for a=1:numel(angles)

    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    Ttiny(:,a)=transmittanceTiny2DCollimated(filter,angles(a),wavelengths,polarization,accuracy);
    

    Tinf(:,a)=transmittanceInfinite(filter,angles(a),wavelengths,polarization);
    
    % -- Ray model --

    %  1. Calculate equivalent monolayer parameters
    normalized_fwhm=fwhm(wavelengths,Tinf(:,1))/targetcwl;  %Normalized FWHM
    R=fwhm2reflectance(normalized_fwhm); % Corresponding mirror reflectance
    
    % 2. Calculate the transmittance
    [T, Tcoh]=transmittanceTinyRayEquivalent(nair,neff,nsub,R,filterwidth,targetcwl,wavelengths,angles(a)+eps,polarization,'accuracy',accuracy,'fastapproximation',false);;
    Tray(:,a)=T;
    Traycoh(:,a)=Tcoh;
    Tray_analytical(:,a)=transmittanceTinyRayEquivalent(nair,neff,nsub,R,filterwidth,targetcwl,wavelengths,angles(a)+eps,polarization,'accuracy',accuracy);   
        
end


%% Plot transmittance
maxnorm=@(x)x/max(x);
maxnorm=@(x)x;
fig=figure(1);clf;  hold on;
fig.Position=[320 168 1273 670];
for a=1:numel(angles)
    subplot(2,3,a); hold on;
    hinf=plot(wavelengths,maxnorm(Tinf(:,a)),':','color','k','linewidth',1.5)
%    htiny(a)=plot(wavelengths,maxnorm(Ttiny(:,a)),'color','k','linewidth',2)
    hray(a)=plot(wavelengths,maxnorm(Tray(:,a)),'-','color',[1 0.2 0.2],'linewidth',2)
    hraycoh(a)=plot(wavelengths,maxnorm(Traycoh(:,a)),'-','color','b','linewidth',2)
    %hrayanalytical(a)=plot(wavelengths,maxnorm(Tray_analytical(:,a)),'--','color',[1 0.2 0.2],'linewidth',2)   
    
    ylabel('Transmittance')
    xlabel('Wavelength (µm)')
    title([num2str(angles(a)) ' deg'])
    box on
end
%legend([htiny(1) hray(1) hrayanalytical(1)   hinf(1)],'Tiny Wave model','Tiny Ray Model','Analytical Ray Model','Infinite filter','location','best')










%% Peak transmittance 
% Demonstrate how the peak transmittance approximately follows the trend
% predicted by considering only ray optics (truncated interference).
% The deviation for the wave model at normal incidence is due to
% diffraction.

figure(2);clf; hold on;


% Calculate relative change in peak transmittance
hray=plot(angles,max(Tray)/max(Tray(:,1)),'linewidth',1.5)
hwave=plot(angles,max(Ttiny)/max(Tray(:,1)),'linewidth',1.5)

% Analytical Expression for peak transmittance

peakAnalytical=peakTransmittanceRay(targetcwl,R,neff,filterwidth,angles)
hanal=plot(angles,peakAnalytical,'k--','linewidth',2)


legend([hanal(1) hray(1) hwave(1)],'Truncated interference','Ray Model','Wave Model')
xlabel('Incidence Angle')
ylabel('Relative transmittance')












% 
% 
% 
%% Vandersluis
clear Iincoh;
for a=1:numel(angles)
    theta=angles(a)+eps; % incidence angle
    thickness=targetcwl/(2*neff);
    
    delta=4*pi*thickness./wavelengths*neff*cosd(theta/neff);
    
    phi=delta;
    
    K=(0.5*filterwidth)/(thickness*tand(theta))
    IA=1./(1+4*R./(1-R)^2    * sin(delta/2).^2);    
    
    %Iincoh(:,a) = IA.*(1+R^2
    % Er scheelt ong iets met deze formule
    Icoh(:,a) = IA.*(1+IA*R^2 ./(K^2*(1-R)^2) .* (1+R.^(2*K)-2*R.^K*cos(delta)) -2*IA.*R/(K*(1-R)^2).*(cos(delta)-R^K*cos((K+1)*delta) -R));
end

figure(2);clf; hold on;
plot(wavelengths,Icoh/max(Icoh(:,1)))
plot(wavelengths,Tray/max(Tray(:,1)),'-','color',[1 0.2 0.2],'linewidth',1)
ylim([0 1])


%% Neuhaus
clear Iincoh;
for a=1:numel(angles)
    theta=angles(a)+eps; % incidence angle
    thickness=targetcwl/(2*neff);
    
    delta=4*pi*thickness./wavelengths*neff*cosd(theta/neff);
    
    
    
    K=filterwidth/(thickness*tand(theta))
    N=K;
    Z=R*exp(1i*delta);
    uneuh(:,a) =1./(1-Z).*(1-1./N.*(Z.*(Z.^(K-N)-Z.^K))./(1-Z));
    Ineuh(:,a) = conj(uneuh(:,a)).*uneuh(:,a);
    
    %Ineuh(:,a) = ((1-R.^K).^2+4*R.^K*sin(K*delta/2).^2) ./ ((1-R).^2+4*R*sin(delta/2).^2);
    
    Ineuh_approx(:,a)=1./((1-R).^2+4*R*sin(delta/2).^2).* (1-(2*R*(1-R)-4*R*sin(delta/2).^2)./(N*((1-R)^2+4*R*sin(delta/2).^2))  +R.^2./(N^2*((1-R)^2+4*R*sin(delta/2).^2)));
end

figure(3);clf; hold on;
plot(wavelengths,Ineuh/max(Ineuh(:,1)),'k');
plot(wavelengths,Traycoh/max(Traycoh(:,1)),'-','color',[1 0.2 0.2],'linewidth',1)
%plot(wavelengths,Ineuh_approx/max(Ineuh_approx(:,1)));
plot(wavelengths,Ttiny/max(Tinf(:,1)),'-','color','b','linewidth',1)
xlim([0.65 0.85])

