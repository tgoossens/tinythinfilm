function [T,Phi_t,Phi_in] =tinytransmittance2d_collimated(filter,incidence_angle,wavelengths,polarization,accuracy)
%  TINYTRANSMITTANCE  Simulate tiny filter transmittance
%   [T] = TINYTRANSMITTANCE2D_collimated(filter,angledeg,wavelengths,polarization,accuracy);
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
%    - Phi_T (Wx1):  Incident flux [W]
%    
%    
%   This function assumes that filter and pixel size are equal.    
%    This function is equivalent to calling TINYTRANSMITTANCE_CORE with wavepacket2d_collimated and a fullwidth pixel kernel.
%    
%    
%  See also TINYFILTER, TINYTRANSMITTANCE_CORE
%    
%  Copyright Thomas Goossens
%  http://github.com/tgoossens

incident_wavepacket = wavepacket2d_collimated(incidence_angle,filter.stack.refractiveindex(1),filter.width);
[T,Phi_t,Phi_in] = tinytransmittance_core(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixel_fullwidth(filter.width));
end