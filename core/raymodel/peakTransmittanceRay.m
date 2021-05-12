function Tpeak = peakTransmittanceRay(central_wavelength,mirror_reflectance,effective_index,filterwidth,angle_deg)
%PEAKTRANSMITTANCERAY Summary of this function goes here
%   Detailed explanation goes here

R= mirror_reflectance;
angle_refracted = asind(sind(angle_deg)/effective_index);
M = min(1e7,effective_index*filterwidth./(central_wavelength*tand(angle_refracted)));
Tpeak = 1+ (1-R.^M).*(3-mirror_reflectance.^M)./(log(R.^(2*M)));
end

