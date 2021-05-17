function Tpeak = peakTransmittanceRay(central_wavelength,mirror_reflectance,effective_index,filterwidth,angle_deg)
%peakTransmittanceRay Calculate the relative drop in peak transmittance for
%a given equivalent tiny monolayer
%
% INPUTS
% ----------
%
% central_wavelength - Central wavelength of the equivalent monolayer
% mirror_reflectance - Mirror reflectances of the equivalent monolayer
% effective_index  - Effective refractive index (cavity of the monolayer)
% filterwidth  - Width of the filter (same units as wavelength)
% angle_deg    - Incidenge angle in degres
%
% OUTPUTS
%----------
% Tpeak - Relative drop in peak amplitude
%    
% The peak amplitude is relative with respect to the infinitely large
% filter (i.e. filterwidth=inf);
%


% Angle of refraction in the monolayer
angle_refracted = asind(sind(angle_deg)/effective_index);

% Maximal number of interfering rays (clipped for large numbers)
M = min(1e7,effective_index*filterwidth./(central_wavelength*tand(angle_refracted))); 

% Peak value formula as derived in the article
Tpeak = 1+ (1-mirror_reflectance.^M).*(3-mirror_reflectance.^M)./(log(mirror_reflectance.^(2*M)));
end

