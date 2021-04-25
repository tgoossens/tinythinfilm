function [T,Phi_t,Phi_in] =transmittanceTiny2DCollimated(filter,incidence_angle_deg,wavelengths,polarization,accuracy)
%  transmittanceTiny2DCollimated  Simulate tiny filter transmittance
%   [T] = transmittanceTiny2DCollimated(filter,angledeg,wavelengths,polarization,accuracy);
%    
%   All inputs are equal o

%   Inputs
%    filter - Struct containing the tiny filter design (See also TINYFILTER)
%    incidence_angle_deg -  Incidence angle  of the collimated light in degrees
%    wavelengths (Wx1) - Wavelengths (same units as filter.width of filter)
%    polarization - ('s' or 'p' or 'unpolarized')    
%    accuracy - 2^floor(accuracy) subdivision of the spatial frequency domain.
%    pixelkernel - Encodes the width of the pixel and which spatial frequencies are sampled 
%
%   Outputs
%    T (Wx1) -  Transmittance of the filter
%    Phi_T (Wx1) -  Transmitted flux [W]
%    Phi_in (Wx1) -  Incident flux [W]
%     
%    
%    This function is equivalent to calling transmittanceTiny2D with
%    wavepacket2d_collimated and a fullwidth pixel kernel
%    (filterwidth=pixelwidth).
%    This is to facilitate a direct function for this very common case.
%    
%    
%  See also tinyfilter, transmittanceTiny2D
%    
%  Copyright Thomas Goossens
%  http://github.com/tgoossens

incident_wavepacket = wavepacket2d_collimated(incidence_angle_deg,filter.stack.refractiveindex(1),filter.width);
[T,Phi_t,Phi_in] = transmittanceTiny2D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixel_fullwidth(filter.width));
end