function eta = admittance(n,wl,nu,polarization) 
% ADMITTANCE Calculate the admittance for each layer for each spatial frequency and wavelength
%   [eta] = ADMITTANCE(n,wl,v,polarization)
%    
%   Inputs
%    - n (Nx1): Stack of refractive indices 
%    - wl (Wx1): Wavelengths    
%    - nu (Vx1): Spatial frequencies nu = k_i cos theta_i, with k_i wavenumber in layer i
%    - polarization ('s' or 'p')    
%    
%   Outputs
%    - @(j)eta : an anonymous function that gives the admittance for a given layer j.  eta(j) returns a matrix of size VxW                 
%
%    
%    
%  The equations in this function were taken from
%  Amra, C., Lequime, M., & Zerrad, M. (2021). Electromagnetic Optics of Thin-Film Coatings: Light Scattering, Giant Field Enhancement, and Planar Microcavities. Cambridge University Press.
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm




% Admittance of free space
    eta0=2.6544E-3; 

    % Correcft dimensions for elementwise multiplcation
    %lambda = reshape(wl,[1 numel(wl)]);
    lambda = reshape(wl,[1 1 numel(wl)]);
    %    nu = reshape(nu,[numel(nu) 1]);

    k = @(j) 2*pi./(lambda).*n(j);
    alpha = @(j) sqrt(k(j).^2-(2*pi*nu).^2);

    
    
    if(polarization == 's')
        eta = @(j) n(j).*eta0.*alpha(j)./k(j);%s polarization
    end

    if(polarization=='p')
        eta = @(j) n(j).*eta0.*k(j)./(alpha(j));%p polarization
    end


end