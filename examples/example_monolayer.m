clear; close all;



%% Create filter


thickness = 10 ; %µm
n = 1.5; % SiO2
nair=1;
nsub=3.56 % silicon
width = 5; % µm
filter=tinyfilter(nair,nsub,n,thickness,width)


%% Choose simulation options

polarisation = 's';

accuracy = 9;
wavelengths=linspace(0.6,0.9,1000); % µm
angles = [0 5 10 15 20]

for a=1:numel(angles)
    

    %% Simulate
    disp(['Simulate tiny filter: ' num2str(angles(a)) ' deg']);
    
    [Ttiny,flux_t]=tinytransmittance(filter,angles(a),wavelengths,polarisation,accuracy);
    T(:,a)=Ttiny;
    flux(:,a)=flux_t;


end


%% Plot transmittance

cmap = hot;
s=size(cmap,1);
color{1}=cmap(1,:);
color{2}=cmap(round(0.45*s),:);
color{3}=cmap(round(0.5*s),:);
color{4}=cmap(round(0.6*s),:);
color{5}=cmap(round(0.66*s),:)


figure(1);clf; 
subplot(211); hold on;
for a=1:numel(angles)
    plot(wavelengths,T(:,a),'color',color{a},'linewidth',2)
end

ylim([0 1])
ylabel('Transmittance')
xlabel('Wavelength (µm)')

subplot(212)
fluxdrop=max(flux,[],1)/max(flux(:,1));
plot(angles,fluxdrop)
ylim([0 1])
ylabel('Relative flux')
xlabel('Incidence angle (deg)')























