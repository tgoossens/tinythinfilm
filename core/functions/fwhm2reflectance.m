function R = fwhm2reflectance(normalized_fwhm)
%FWHM2REFLECTANCE Transforms normalized fwhm to equivalent mirror
%reflectance for the equivalent monolayer model
%
% INPUT
%     normalized_fwhm  =   (FWHM)/(central wavelength)
% OUTPUT
%    R = Mirror reflectance (between zero and 1)
R = 2-cos(pi*normalized_fwhm)-sqrt(3-4*cos(pi*normalized_fwhm)+cos(normalized_fwhm*pi)^2);

end

