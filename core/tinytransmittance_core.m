function [T,Phi_t,Phi_in] = tinytransmittance_core(filter,angledeg,wavelengths,polarization,accuracy,pixelkernel);
%  TINYTRANSMITTANCE  Simulate tiny filter transmittance
%   [T] = TINYTRANSMITTANCE(filter,angledeg,wavelengths,polarization,accuracy);
%    
%   Inputs
%    - filter : Struct containing the tiny filter design (See also TINYFILTER)
%    - angledeg:  Incidence angle in degrees
%    - wavelengths (Wx1): Wavelengths (same units as filter.width of filter)
%    - polarization ('s' or 'p')    
%    - accuracy: 2^floor(accuracy) subdivision of the spatial frequency domain.
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - Phi_T (Wx1):  Transmitted flux [W]
%    - Phi_T (Wx1):  Incident flux [W]
%    
%    
%  See also TINYFILTER    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens
    
    
wl=wavelengths;
anglerad=deg2rad(angledeg); 
width=filter.width;


% Spatial frequency integration domain
nu = linspace(-1/wl(1), 1/wl(1),2^floor(accuracy))';



%% Definitions and helper ufcntions

% Wavenumber
k = @(n) 2*pi./(wl)*n; 

% Fourier transform of the pixel kernel (so we don't recompute it for each wavelength)
fftpix=fft(pixelkernel(nu));

%%  Calculate admittances

% Complex surface admittance of filter stack
% We will only use the transmission coefficient here
[Y0,r,t] = surfaceadmittance(filter.n,filter.h,wl,nu,polarization);

% Admittances of each layer
eta = admittance(filter.n,wl,nu,polarization);
eta_in=eta(1);
eta_sub=eta(numel(filter.n));

for j=1:numel(wl)
    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wwave
    Ain(:,j) = width*sinca(pi*width*(nu-filter.n(1)*sin(anglerad)/wl(j))); 
    
    % Useful integration domain;. This conditions corresponds to ignore incidence angles larger than 90 degres.
    domain = abs(nu).*wl(j) <=1;
    
    % Transmitted wave
    At(:,j)=domain.*t(:,j).*Ain(:,j);


    
    %%%%%%%%%% FLUXES  %%%%%%%%%%%%
    %Incident flux
    temp = real(eta_in(:,j).*Ain(:,j).*fftshift(ifft(fft(conj(Ain(:,j))).*fftpix)));
    Phi_in(j)=trapz(nu,temp);

    
    % Transmitted flux
    temp=real(eta_sub(:,j).*At(:,j).*fftshift(ifft(fft(conj(At(:,j))).*fftpix)));
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

