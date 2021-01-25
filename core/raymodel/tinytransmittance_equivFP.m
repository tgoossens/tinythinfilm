function [T] = tinytransmittance_equivFP(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,polarization,accuracy) 
% TINYTRANSMITTANCE_EQUIVFP
% Simulate tiny transmittance of a monolayer using a ray based model    
%
% [T] = TINYTRANSMITTANCE_EQUIVFP(n0, nlayer,nsubstrate,height,width,wavelengths,angles,polarization,accuracy)
%     
%   Inputs
%    - n0: Refractive index incident medium
%    - nlayer: Refractive index monolayer
%    - nsubstrate:Refractive index substrate
%    - R: Product of reflection coefficients
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
    equivstack = [1 neff 1];
    
% Cosineof refraction angle angle
    costh_n = sqrt(1-sind(angledeg).^2./equivstack.^2);
    
    % Calculate admittances depending on polarization    
    if(polarization=='s')    
        eta0=n0*costh_n(1);        
        eta1=neff*costh_n(2);        
        eta2=nsub*costh_n(3);
    elseif(polarization=='p')   
        eta0=n0/costh_n(1);        
        eta1=neff/costh_n(2);        
        eta2=nsub/costh_n(3);
    end
    
    % Mirror properties
    t01=2*eta0/(eta0+eta1)
    t12=2*eta1/(eta1+eta2)
    t02=2*eta0/(eta0+eta2)


    % Filter dimensions
    height=cwl/(2*neff);
    w=width;
    
    % Sampling of spatial pixel axis
    x=linspace(0,width,2^floor(accuracy));
    dx = abs(x(2)-x(1));
    
    % Calculation of number of reflections (anonymous fucntion)
    th_n = @(th) asind(sind(th)/equivstack(2));
    num=@(x,th)x./(height*tand(th_n(th)));
    M = @(x,th) floor(num(x,th)/2+1/2)
    
    T = zeros(numel(wavelengths),1);
    delta = 2*pi*equivstack(2)*height*costh_n(2)./wavelengths';
    for ix=1:numel(x)
        m = M(x(ix),angledeg);  % number of rays
        formula=(R^(2*m)-2*R^(m)*cos(2*m*delta)+1)./(R^(2)-2*R*cos(2*delta)+1);
        T = T+dx*formula;
    end
    
    T = (1-R)^2*T/width;
end
