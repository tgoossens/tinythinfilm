function [T,Phi_t,Phi_in] =transmittanceTiny2DCollimated(filter,incidence_angle,wavelengths,polarization,accuracy)
%  transmittanceTiny2DCollimated  Simulate tiny filter transmittance
%   [T] = transmittanceTiny2DCollimated(filter,angledeg,wavelengths,polarization,accuracy);
%    
%   Inputs
%    - filter : Struct containing the tiny filter design (See also TINYFILTER)
%    - angledeg:  Incidence angle in degrees
%    - wavelengths (Wx1): Wavelengths (same units as filter.width of filter)
%    - polarization ('s' or 'p' or 'unpolarized')    
%    - accuracy: 2^floor(accuracy) subdivision of the spatial frequency domain.
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - Phi_T (Wx1):  Transmitted flux [W]
%    - Phi_in (Wx1):  Incident flux [W]
%    
%    
%   This function assumes that filter and pixel size are equal.    
%    This function is equivalent to calling TINYTRANSMITTANCE_CORE with wavepacket2d_collimated and a fullwidth pixel kernel.
%    
%    
%  See also tinyfilter, transmittanceTiny2D
%    
%  Copyright Thomas Goossens
%  http://github.com/tgoossens

incident_wavepacket = wavepacket2d_collimated(incidence_angle,filter.stack.refractiveindex(1),filter.width);
[T,Phi_t,Phi_in] = transmittanceTiny2D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixel_fullwidth(filter.width));
end