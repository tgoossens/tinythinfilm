function [T,Phi_t,Phi_in] = classictransmittance_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,polarization,accuracy);
%  CLASSICTRANSMITTANCE_MONO  Simulate infinitely wide transmittance of an equivalent monolayer with a given bandwidth, size and effective index.
%   [T] = CLASSICTRANSMITTANCE_MONO(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,polarization,accuracy);
%    
%   Inputs
%    - central_wavelength: Central wavelength of the filter
%    - normalized_fwhm : FWHM as a percentage of central wavelength (e.g.  0.01), called  Lambda_Infty in the associated paper. This is the FWHM at normal incidence.
%    - effective_index: The effective refractive index which determines angular sensitivity
%    - width: Spatial size of the filter    
%    - nsub: Substrate refractive index (for proper scaling of the peak transmittance)
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
%  See also TINYTRANSMITTANCE_MONO  CLASSICTRANSMITTANCE    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens

    
    
    
    if(or(polarization=='unpolarized',polarization=='unpolarised'))
        [T_s,Phi_t_s,Phi_in_s] = classictransmittance_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,'s',accuracy);
        [T_p,Phi_t_p,Phi_in_p] = classictransmittance_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,'p',accuracy);
        T =  0.5*(T_s+T_p);
        Phi_t =  0.5*(Phi_t_s+Phi_t_p);
        Phi_in =  0.5*(Phi_t_s+Phi_t_p);
        
        return;
    end
    
%% Renaming of some variables
wl=wavelengths;
anglerad=deg2rad(angledeg); 
neff=effective_index;
cwl=central_wavelength;



%% Filter size
% The cavity thickness is chosen to be a halve-wave plate such that constructive interference occurs at the desired wavelength.
h = cwl/(2*neff);

%% Conversion functions
d=@(alpha) 0.5*pi*alpha; %delta
fwhm2r=@(alpha)-sqrt(cos(2*d(alpha)).^2-4*cos(2*d(alpha))+3)-cos(2*d(alpha))+2    ;


% Convert the normalized_fwhm to equivalent mirror reflectance
R = fwhm2r(normalized_fwhm);

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


%% Transmittance from incident medium to substrate
ts=2*eta_in./(eta_sub+eta_in); 
Ts=conj(ts).*ts;

%% Transmittance of the infinitely wide filter
% Transmitted flux
Phi_t = real(eta_sub)*Ts*(1-R)^2/2 * 1./(1-2*R*cos(2*delta)+R^2);

% Incident flux
Phi_in = real(eta_in)/2;

% Transmittance
T=Phi_t./Phi_in;


function f = sinca(x)
% Modified sinc function because matlab sinc function already includes the factor pi.
% This makes notation consistent with definitions in the publications.    
    f=sinc(x/pi); 
end
end

