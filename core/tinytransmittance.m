function [T,Phi_t,Phi_in] = tinytransmittance(filter,angledeg,wavelengths,polarization,accuracy);
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
%   This function assumes that filter and pixel size are equal.    
%    This function is equivalent to calling TINYTRANSMITTANCE_CORE with a fullwidth pixel kernel.
%    
%    
%  See also TINYFILTER, TINYTRANSMITTANCE_CORE
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens
    
    [T,Phi_t,Phi_in] = tinytransmittance_core(filter,angledeg,wavelengths,polarization,accuracy,pixel_fullwidth(filter.width));
end