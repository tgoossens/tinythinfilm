function [T] = tinytransmittance_ray(n0,nlayer,nsubstrate,height,width,wavelengths,angledeg,polarization,accuracy)
% TINYTRANSMITTANCE_RAY
% Simulate tiny transmittance of a monolayer using a ray based model    
%
% [T] = TINYTRANSMITTANCE_RAY(n0, nlayer,nsubstrate,height,width,wavelengths,angles,polarization,accuracy)
%    
%   Inputs
%    - n0: Refractive index incident medium
%    - nlayer: Refractive index monolayer
%    - nsubstrate:Refractive index substrate
%    - height: Thickness of the film
%    - width: Width of the film (same units as height)
%    - wavelengths (Wx1): Wavelengths (same units as height)
%    - angledeg:  Incidence angle in degrees
%    - polarization ('s' or 'p')    
%    - accuracy: Subdivide integration domain in 2^floor(accuracy) points
%   Outputs
%    - T (Wx1):  Transmittance of tiny monolayer
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens
%    

% Full thin film stack
    stack = [n0 nlayer nsubstrate];
    
% Cosineof refraction angle angle
    costh_n = sqrt(1-sind(angledeg).^2./stack.^2);
    
    % Calculate admittances depending on polarization    
    if(polarization=='s')    
        eta0=stack(1)*costh_n(1);        
        eta1=stack(2)*costh_n(2);
        eta2=stack(3)*costh_n(3);
    elseif(polarization=='p')   
        eta0=stack(1)/costh_n(1);        
        eta1=stack(2)/costh_n(2);
        eta2=stack(3)/costh_n(3);
    end

    % Mirror properties
    t01=2*eta0/(eta0+eta1)
    t12=2*eta1/(eta1+eta2)
    r10=(eta1-eta0)/(eta0+eta1);
    r12=(eta1-eta2)/(eta1+eta2);
    r=r10*r12;
    
    % Filter dimensions
    h=height;
    w=width;
    
    % Sampling of spatial pixel axis
    x=linspace(0,width,2^floor(accuracy));
    dx = abs(x(2)-x(1));
    
    % Calculation of number of reflections (anonymous fucntion)
    th_n = @(th) asind(sind(th)/stack(2));
    num=@(x,th)x./(h*tand(th_n(th)));
    M = @(x,th) floor(num(x,th)/2+1/2)
    
    T = zeros(numel(wavelengths),1);
    delta = 2*pi*stack(2)*height*costh_n(2)./wavelengths';
    for ix=1:numel(x)
        m = M(x(ix),angledeg);  % number of rays
        formula=(r^(2*m)-2*r^(m)*cos(2*m*delta)+1)./(r^(2)-2*r*cos(2*delta)+1);
        T = T+dx*formula;
    end
    
    T = ((t01*t12)^2*(eta2/eta0) *T)/width;
end
