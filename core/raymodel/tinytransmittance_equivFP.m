function [T] = tinytransmittance_equivFP(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,polarization,accuracy) 
% TINYTRANSMITTANCE_EQUIVFP
% Simulate tiny transmittance of a Fabry-Pérot using a ray based model    
%
% [T] = TINYTRANSMITTANCE_EQUIVFP(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,polarization,accuracy) 
%     
%   Inputs
%    - n0: Refractive index incident medium
%    - neff: effective refractive index of the cavity
%    - nsubstrate:Refractive index substrate
%    - R: Product of reflection coefficients
%    - width: Width of the film
%    - cwl: Central wavelength of the filter (same units as width)
%    - wavelengths (Wx1): Wavelengths (same units as height)
%    - angledeg:  Incidence angle in degrees
%    - polarization ('s' or 'p')    
%    - accuracy: Subdivide integration domain in 2^floor(accuracy) points
%   Outputs
%    - T (Wx1):  Ray-model estimation of the transmittance of a tiny Fabry-Pérot filter 
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm
%    

% Full thin film stack
    equivstack = [1 neff 1];
    
    % Cosineof refraction angle
    costh_n = sqrt(1-sind(angledeg).^2./equivstack.^2);
    
    % Calculate characteristic admittances depending on polarization    
    if(polarization=='s')    
        eta0=n0*costh_n(1);        
        eta1=neff*costh_n(2);        
        eta2=nsub*costh_n(3);
    elseif(polarization=='p')   
        eta0=n0/costh_n(1);        
        eta1=neff/costh_n(2);        
        eta2=nsub/costh_n(3);
    end
    

    % Filter dimensions
    height=cwl/(2*neff);
    w=width;
    
    % Sampling of spatial pixel axis
    x=linspace(0,width,2^floor(accuracy));
    dx = abs(x(2)-x(1));
    
    % Calculation of number of reflections (anonymous fucntion)
    th_n = @(th) asind(sind(th)/equivstack(2));
    num=@(x,th)x./(height*tand(th_n(th)));
    N = @(x,th) floor(num(x,th)/2+1/2)
    
    T = zeros(numel(wavelengths),1);
    delta = 2*pi*equivstack(2)*height*costh_n(2)./wavelengths';
    for ix=1:numel(x)
        n = N(x(ix),angledeg);  % number of rays
        formula=(R^(2*n)-2*R^(n)*cos(2*n*delta)+1)./(R^(2)-2*R*cos(2*delta)+1);
        T = T+dx*formula;
    end
    
    Tsub = 1-((1-nsub)/(1+nsub))^2;
    T = Tsub*(1-R)^2*T/width;
end
