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
    
    
wl=reshape(wavelengths,[1 1 numel(wavelengths)]);
anglerad=deg2rad(angledeg); 
width=filter.width;


% Spatial frequency integration domain
% Exact evaluation at -1/wl(1) would result in a divison by zero for the p-polarized calculation
nu = linspace(-0.99/wl(1), 0.99/wl(1),2^floor(accuracy))';
nu = reshape(nu,[numel(nu) 1 1]);
%nu = linspace(-10/wl(1), 10/wl(1),2^floor(accuracy))';



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
[Y0,r,t] = surfaceadmittance(filter.n,filter.h,wl,nu,polarization);

% Admittances of each layer
eta = admittance(filter.n,wl,nu,polarization);
eta_in=eta(1);
eta_sub=eta(numel(filter.n));

for j=1:numel(wl)
    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wwave
    Ain(:,1,j) = width*sinca(pi*width*(nu-filter.n(1)*sin(anglerad)/wl(j))); 
    

    
    % Useful integration domain;. This conditions corresponds to ignore incidence angles larger than 90 degres.
    domain = abs(nu).*wl(j) <=1;
    
    % Transmitted wave
    At(:,1,j)=domain.*t(:,j).*Ain(:,j);


    
    %%%%%%%%%% FLUXES  %%%%%%%%%%%%
    %Incident flux (explicit result)
    nu_angle=filter.n(1)*sin(anglerad)/wl(j);
    eta_in = admittance(filter.n,wl(j),nu_angle,polarization);
    Phi_in(j)=real(eta_in(1))/2 * filter.width;
    
     
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

