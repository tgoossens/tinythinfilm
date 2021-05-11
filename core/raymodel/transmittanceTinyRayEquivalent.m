function [T] = transmittanceTinyRayEquivalent(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,polarization,accuracy,flag_fastapproximation) 
% transmittanceTinyRayEquivalent
% Simulate tiny transmittance of a Fabry-Pérot using a ray based model    
%
% [T] = transmittanceTinyRayEquivalent(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,polarization,accuracy) 
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
%    - flag_fastapproximation: Calculate the transmittance using an
%    analytical approxmiation that evaluates faster. By default false.
%   Outputs
%    - T (Wx1):  Ray-model estimation of the transmittance of a tiny Fabry-Pérot filter 
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm
%    


 if ~exist('flag_fastapproximation','var')
     % Default value:
      flag_fastapproximation=false;
 end


if(or(polarization=='unpolarized',polarization=='unpolarised'))
    [T_s] = transmittanceTinyRayEquivalent(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,'s',accuracy,flag_fastapproximation) ;
    [T_p] = transmittanceTinyRayEquivalent(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,'p',accuracy,flag_fastapproximation) ;
    T =  0.5*(T_s+T_p);
    return;
end


% Full thin film stack
    equivstack = [1 neff nsub];
    
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
        
    % Transmittance between incident medium and substrate (in absence of
    % monolayer)
    Tsub = 1-((eta0-eta2)/(eta0+eta2))^2;         
    
    % Filter dimensions
    height=cwl/(2*neff);
    w=width;
    
    % Sampling of spatial pixel axis
    x=linspace(0,width,2^floor(accuracy));
    dx = abs(x(2)-x(1));
    
    % Calculation of number of reflections (anonymous fucntion)
    th_n = @(th) asind(sind(th)/equivstack(2));
    num=@(x,th)x./(height*tand(th_n(th)));
    N = @(x,th) floor(num(x,th)/2+1/2);
    
    T = zeros(numel(wavelengths),1);
    delta =(2*pi*equivstack(2)*height*costh_n(2)./wavelengths') - pi; % (because of continous approximation, see supplementary material)
    
    if(flag_fastapproximation)
        % Analytical approximaation
        M = floor(min(width*neff/(cwl*tand(angledeg/neff)),1e7));
        T=  Tsub.*(1-R).^2 .* (1+ (R.^(2*M)-1)./log(R.^(2*M))-2*(log(R).*(R.^M .*cos(2*M*delta)-1)+2*delta.*R.^M.*sin(2*M.*delta))./(4*M.*delta.^2+M.*log(R).^2))./(1-2*R*cos(2*delta)+R.^2);
    
    else
        for ix=1:numel(x)
            n = min(N(x(ix),angledeg),1e7);  % Number of interfering rays is limited to avoid division by zero when n=Inf
            formula=(R^(2*n)-2*R^(n)*cos(2*n*delta)+1)./(R^(2)-2*R*cos(2*delta)+1);
            T = T+dx*formula;
        end
        T = Tsub*(1-R)^2*T/width;
    end


    
    

 
end
