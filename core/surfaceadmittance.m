function [Y0,r,t] = surfaceadmittance(n,h,wl,nu,polarization) 
% SURFACEADMITTANCE  Calculate the complex surface admittance and the corresponding transmission and reflection coefficients
%   [Y,r,t] = SURFACEADMITTANCE(n,h,wl,v,polarization)
%    
%   Inputs
%    - n (Nx1): Stack of refractive indices 
%    - h (Nx1): Stack of layer thicknesses
%    - wl (Wx1): Wavelengths    
%    - nu (Vx1): Spatial frequencies nu = k_i cos theta_i, with k_i wavenumber in layer i
%    - polarization ('s' or 'p')    
%   Outputs
%    - Y0 (VxW):  Complex surface admittance for each wavelength and spatial frequency
%    - r (VxW):  Reflection coefficient for each wavelength and spatial frequency
%    - t (VxW):  Transmission coefficient for each wavelength and spatial frequency
%
%    
%    
%  The equations in this function were taken from
%  Lequime, M., and Amra, C. De l’Optique électromagnétique à l’Interférométrie-Concepts et illustrations: Concepts et illustrations. EDP Sciences, 2013
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens




% Admittance of free space
    eta0=2.6544E-3; 

    % Correcft dimensions for elementwise multiplcation
    lambda = reshape(wl,[1 numel(wl)]);
    nu = reshape(nu,[numel(nu) 1]);

    k = @(j) 2*pi./(lambda).*n(j);
    alpha = @(j) sqrt(k(j).^2-(2*pi*nu).^2);

    if(polarization == 's')
        eta = @(j) n(j).*eta0.*alpha(j)./k(j);%s polarization
    end

    if(polarization=='p')
        eta = @(j) n(j).*eta0.*k(j)./(alpha(j));%p polarization
    end

    % Phase thickness of each layer
    delta = @(j) alpha(j).*h(j);

    N = numel(n);

    %% initiate recursion
    Y = eta(N);

    tfactor=1;

    j=2;
    for j=(N-1):-1:2
        tfactor = tfactor.*(cos(delta(j))-1i*Y./eta(j).* ...
                            sin(delta(j))); 
        
        teller=(Y-1i*eta(j).*tan(delta(j)));    
        noemer=(1-1i*(Y./eta(j)).*tan(delta(j)));
        Y =  teller./noemer;
    end


    r = (eta(1)-Y)./(eta(1)+Y);
    t = (1+r)./tfactor;

    % Pageina 108 in amra

    % Y0 is the complex surface admittance
    Y0=Y;

end