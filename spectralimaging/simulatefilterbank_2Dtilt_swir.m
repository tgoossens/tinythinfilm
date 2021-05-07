clear;

%% Define filters
nbFilters = 60;
nbWavelengths=2^9;
wavelengths=linspace(1000,1700,nbWavelengths);

cwl = linspace(1100,1650,nbFilters);
nFWHM=10/1100
width=15000; %nm


filterwidth=width;

neff=1.7;
nsub=3.67;%silicon
nsub=1;
nair=1;



%% Simulation parameters
polarization = 's'
accuracy=9

angledeg = [0.1 10 15];
fnumber=6;
coneangle_deg = atand(1./(2*fnumber));
alpha_deg=0; % azimuth

realwidth=sqrt(0.42)*5.5
pixelkernel=pixel_partialwidth(-realwidth/2,realwidth/2);
pixelkernel=pixel_fullwidth(filterwidth);

%%
tic
for a=1:numel(angledeg)
        for c=1:nbFilters
         %disp(['CRA ' num2str(angledeg(a)) ' deg - Filter ' num2str(f)]),
         filter=tinyfilterCreateEquivalent(cwl(c),nFWHM,neff,width,nair,nsub);
         Tinf(:,c,a)=transmittanceInfinite(filter,angledeg(a),wavelengths,polarization);

        end
end
toc

%% Simulation incidence angles
tic
for a=1:numel(angledeg)

        for c=1:nbFilters
         filter=tinyfilterCreateEquivalent(cwl(c),nFWHM,neff,width,nair,nsub);
         wavepacket2d = wavepacket2DCollimated(angledeg(a),nair,filterwidth);
          Twave(:,c,a)=transmittanceTiny2DCollimated(filter,angledeg(a),wavelengths,polarization,accuracy);
          T(:,c,a)=transmittanceTinyRayEquivalent(nair,neff,nsub,1-pi*nFWHM,width,cwl(c),wavelengths,angledeg(a),polarization,accuracy);
        end

end
toc


%% Reflectance of sample
clear R;
cyclespernm=0.01;
R = 1+sin(2*pi*wavelengths*cyclespernm)';
%% SImulate flat fielding
% White reference

white = ones(size(R));

%  Calculate sensor output §with flat fielding

for a=1:numel(angledeg)
   %Filter response matrix for a specific angle
   B = T(:,:,a);    
   Binf = Tinf(:,:,a);
   
   % Flat fielding
   d(:,a)= (B'*R)  ./ (B'*white);
   dinf(:,a)= (Binf'*R)./(Binf'*white);
end



% Plot
fig=figure(1);clf;hold on;
fig.Position=[610 87 988 721]
for a=1:numel(angledeg)
    subplot(numel(angledeg),1,a); hold on;
    cwlcorr=cwl.*(1-deg2rad(coneangle_deg).^2./(4*neff^2)-deg2rad(angledeg(a)).^2/(2*neff^2));
    plot(cwlcorr,d(:,a),'k.-')
%    plot(cwlcorr,dwp(:,a),'b.-')
    plot(cwlcorr,dinf(:,a),'r.-')
    plot(wavelengths,R,'k:')
    
    title(['f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg'])
    xlabel('Wavelength (nm)')
    ylabel('Transmittance(%)')
    %plot(wavelengths,T(:,:,a))
    %plot(wavelengths,Tinf(:,:,a),'k:')
end

legend('tiny ray','inf','reference','location','best')
%%
%saveas(gcf,['./fig/swir_collimated_width' num2str(width/1000) '-nfwhm' num2str(nFWHM*100) 'percent-brownglass.png'])




%% Plot filter respones

figure(2);clf;hold on;
for a=1:numel(angledeg)
   subplot(numel(angledeg),1,a);
   plot(wavelengths,T(:,:,a))
       title(['f/' num2str(fnumber) ' - CRA ' num2str(angledeg(a)) ' deg'])
       ylim([0 1])
    xlabel('Wavelength (nm)')
    ylabel('Transmittance(%)')
end




