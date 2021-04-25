function [T,R,t,r] = infinitetransmittance(filter,angledeg,wavelengths,polarization);
%  CLASSICTRANSMITTANCE  Simulate filter transmittance for an infinitely wide filter
%   [T,R,t,r] = CLASSICTRANSMITTANCE(filter,angledeg,wavelengths,polarization);
%    
%   Inputs
%    - filter : Struct containing the tiny filter design (See also TINYFILTER)
%    - angledeg:  Incidence angle in degrees
%    - wavelengths (Wx1): Wavelengths (same units as filter.width of filter)
%    - polarization ('s' or 'p' or 'unpolarized')    
%
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - R (Wx1):  Reflectance of the filter
%    - t (Wx1):  Transmission coefficient
%    - r (Wx1):  Reflection coefficient of the filter
%    
%  See also TINYTRANSMITTANCE, TINYFILTER    
%    
%  The method implemented can be found in
%    Lequime, M., and Amra, C. De l’Optique électromagnétique à
% l’Interférométrie-Concepts et illustrations: Concepts et illustrations. EDP
%  Sciences, 2013
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens



    if(or(polarization=='unpolarized',polarization=='unpolarised'))
        [T_s,R_s,t_s,r_s] = infinitetransmittance(filter,angledeg,wavelengths,'s');
        [T_p,R_p,t_p,r_p] = infinitetransmittance(filter,angledeg,wavelengths,'p');
        T =  0.5*(T_s+T_p);
        R =  0.5*(R_s+R_p);
        t =  0.5*(t_s+t_p);
        r =  0.5*(r_s+r_p);
        
        return;
    end
    
    %%  Calculate admittances

    % Complex surface admittance of filter stack
    % We will only use the transmission coefficient here
    
    for w=1:numel(wavelengths)
        k = ((2*pi)./wavelengths(w)) * filter.stack.refractiveindex(1); 
        sigma = k*sind(angledeg); nu=sigma/(2*pi);
        t=filter.transmission(wavelengths(w),nu,polarization);
        r=NaN;
        % Admittances of each layer
        eta = admittance(filter.stack.refractiveindex,wavelengths(w),nu,polarization);
        eta_in=eta(1);
        eta_sub=eta(numel(filter.stack.refractiveindex));
        
        
        % Transmittance and Reflectance can be calculated form transmission and reflecction coefficients
        % Take real part to remove numerical imaginary residue
        T(w) = real(eta_sub)./real(eta_in) .*  real(conj(t).*t);
        R(w) = NaN; % TODO
        
    end
end

