function [T,Phi_t,Phi_in] = classictransmittance_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,polarization,accuracy);
%  TINYTRANSMITTANCE  Simulate tiny filter transmittance
%   [T] = TINYTRANSMITTANCE_CORE_MONO(filter,angledeg,wavelengths,polarization,accuracy);
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
neff=effective_index;
cwl=central_wavelength;



%% Definitions alpha and k

% Half wave plate
h = cwl/(2*neff);

%% Conversion functions
d=@(alpha) 0.5*pi*alpha; %delta
fwhm2r=@(alpha)-sqrt(cos(2*d(alpha)).^2-4*cos(2*d(alpha))+3)-cos(2*d(alpha))+2    



R = fwhm2r(normalized_fwhm)

%%  Calculate admittances

% Complex surface admittance of filter stack
% We will only use the transmission coefficient here
% Admittances of each layerT

delta = 2*pi*neff*h*sqrt(1-sind(angledeg)^2/neff^2)./wl;

nu_angle=sin(anglerad)/wl(1);
eta = admittance(nsub,wl(1),nu_angle,polarization);
eta_sub=eta(1);
eta = admittance(1,wl(1),nu_angle,polarization);
eta_in=eta(1);


ts=2*eta_in./(eta_sub+eta_in); 
Ts=conj(ts).*ts;
% Transmitted flux
Phi_t = real(eta_sub)*Ts*(1-R)^2/2 * 1./(1-2*R*cos(2*delta)+R^2);
Phi_in = real(eta_in)/2;


% Transmittance
T=Phi_t./Phi_in;





function f = sinca(x)
% Modified sinc function because matlab sinc function already includes the factor pi.
% This makes notation consistent with definitions in the publications.    
    f=sinc(x/pi); 
end
end

