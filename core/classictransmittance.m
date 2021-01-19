function [T,R,t,r] = classictransmittance(filter,angledeg,wavelengths,polarization);
%  CLASSICTRANSMITTANCE  Simulate filter transmittance for an infinitely wide filter
%   [T,R,t,r] = CLASSICTRANSMITTANCE(filter,angledeg,wavelengths,polarization,accuracy);
%    
%   Inputs
%    - filter : Struct containing the tiny filter design (See also TINYFILTER)
%    - angledeg:  Incidence angle in degrees
%    - wavelengths (Wx1): Wavelengths (same units as filter.width of filter)
%    - polarization ('s' or 'p')    
%
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - R (Wx1):  Reflectance of the filter
%    - t (Wx1):  Transmission coefficient
%    - r (Wx1):  Reflection coefficient of the filter
%    
%  See also TINYTRANSMITTANCE, TINYFILTER    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens



    %%  Calculate admittances

    % Complex surface admittance of filter stack
    % We will only use the transmission coefficient here
    
    for w=1:numel(wavelengths)
        k = ((2*pi)./wavelengths(w)) * filter.n(1); 
        sigma = k*sind(angledeg); nu=sigma/(2*pi);
        [Y0,r,t] = surfaceadmittance(filter.n,filter.h,wavelengths(w),nu,polarization);

        % Admittances of each layer
        eta = admittance(filter.n,wavelengths(w),nu,polarization);
        eta_in=eta(1);
        eta_sub=eta(numel(filter.n));
        
        
        % Transmittance and Reflectance can be calculated form transmission and reflecction coefficients
        T(w) = real(eta_sub./eta_in .*  conj(t)*t);
        R(w) = conj(r)*r;
        
    end
end

