function normalized_fwhm = reflectance2fwhm(R)
%FWHM2REFLECTANCE Transforms equivalent mirror transmittance to normalized
%fwhm as used in the context of the equivalent monolayer model
%
% INPUT
%    R  - Mirror reflectance (between zero and 1)
% OUTPUT
%     normalized_fwhm  =   (FWHM)/(central wavelength)

normalized_fwhm = 1/pi * acos(2-R/2-1/(2*R));

end

