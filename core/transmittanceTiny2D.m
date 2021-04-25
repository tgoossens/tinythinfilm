function [T,Phi_t,Phi_in] = transmittanceTiny2D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixelkernel)
% function [T,Phi_t,Phi_in] = transmittanceTiny2D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixelkernel)
%  transmittanceTiny2D  Simulate tiny filter transmittance
%    
%   Inputs
%    filter - Struct containing the tiny filter design (See also TINYFILTER)
%    incident_wavepacket -  A function @(nu,wavelength) of spatial frequency and wavelength that describes the angular spectrum of the field that enters the filter
%                           see directory wavepackets/
%    wavelengths (Wx1) - Wavelengths (same units as filter.width of filter)
%    polarization - ('s' or 'p' or 'unpolarized')    
%    accuracy - 2^floor(accuracy) subdivision of the spatial frequency domain.
%    pixelkernel - Encodes the width of the pixel and which spatial frequencies are sampled 
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - Phi_T (Wx1):  Transmitted flux [W]
%    - Phi_in (Wx1):  Incident flux [W]
%    
%    
%  See also TINYFILTER    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens

    
    if(or(polarization=='unpolarized',polarization=='unpolarised'))
        [T_s,Phi_t_s,Phi_in_s] = transmittanceTiny2D(filter,incident_wavepacket,wavelengths,'s',accuracy,pixelkernel);
        [T_p,Phi_t_p,Phi_in_p] = transmittanceTiny2D(filter,incident_wavepacket,wavelengths,'p',accuracy,pixelkernel);
        T =  0.5*(T_s+T_p);
        Phi_t =  0.5*(Phi_t_s+Phi_t_p);
        Phi_in=  0.5*(Phi_in_s+Phi_in_p);
        return;
    end
    
    
wl=reshape(wavelengths,[1 1 numel(wavelengths)]);
width=filter.width;


% Spatial frequency integration domain
% Exact evaluation at -1/wl(1) would result in a divison by zero for the p-polarized calculation
nu = linspace(-0.99/wl(1), 0.99/wl(1),2^floor(accuracy))';
dnu = (1/(2*width));


nu = reshape(nu,[numel(nu) 1 1]);



%% Definitions and helper ufcntions

% Wavenumber
k = @(n) 2*pi./(wl)*n; 

% Fourier transform of the pixel kernel (so we don't recompute it for each wavelength)
fftpix=fft(pixelkernel(nu));
conv_pix=@(f) conv(f,pixelkernel(nu),'same');
% To be implemented: convolution using FFT


%%  Calculate admittances 

% Complex surface admittance of filter stack
% We will only use the transmission coefficient here
t=filter.transmission(wl,nu,polarization);
% Admittances of each layer
eta = admittance(filter.stack.refractiveindex,wl,nu,polarization);
eta_in=eta(1);
eta_sub=eta(numel(filter.stack.refractiveindex));

for j=1:numel(wl)
    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wwave
    %Ain(:,1,j) = width*sinca(pi*width*(nu-filter.stack.refractiveindex(1)*sin(anglerad)/wl(j))); 
    incident_wavepacket(nu,wl(j));
    Ain(:,1,j) = incident_wavepacket(nu,wl(j));

    
    % Useful integration domain;. This conditions corresponds to ignore incidence angles larger than 90 degres.
    domain = abs(nu).*wl(j) <=1;
    
    % Transmitted wave
    At(:,1,j)=domain.*t(:,j).*Ain(:,j);


    
    %%%%%%%%%% FLUXES  %%%%%%%%%%%%
    %Incident flux (explicit result)
%     nu_angle=filter.stack.refractiveindex(1)*sin(anglerad)/wl(j);
%     eta_in = admittance(filter.stack.refractiveindex(1),wl(j),nu_angle,polarization);
%     Phi_in(j)=real(eta_in(1))/2 * filter.width;
%     
%     
    % Incident flux
    temp=  0.5*real(eta_in(:,1,j).*Ain(:,1,j).*conv_pix(conj(Ain(:,1,j))));
    temp= temp*abs(nu(2)-nu(1)); % discretization convolution integral
    Phi_in(j)=trapz(nu,temp); 
     
    % Transmitted flux
    temp=  0.5*real(eta_sub(:,1,j).*At(:,1,j).*conv_pix(conj(At(:,1,j))));
    temp= temp*abs(nu(2)-nu(1)); % discretization convolution integral
    Phi_t(j)=trapz(nu,temp); 
    
end


% Transmittance
T=Phi_t./Phi_in;




function f = sinca(x)
% Modified sinc function because matlab sinc function already includes the factor pi.
% This makes notation consistent with definitions in the publications.    
    f=sinc(x/pi); 
end
end

