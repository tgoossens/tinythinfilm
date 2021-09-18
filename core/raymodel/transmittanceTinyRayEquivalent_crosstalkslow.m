function [T] = transmittanceTinyRayEquivalent_crosstalk(n0,neff,nsub,R,filterwidths,cwls,wavelengths,angledeg,polarization,accuracy,pixelrange,coherent_flag)
% transmittanceTinyRayEquivalent
% Simulate tiny transmittance of a Fabry-Pérot using a ray based model
%
% [T] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,polarization,accuracy,pixelrange,flag_fastapproximation)
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
%
%  Variable inputs
%    - 'accuracy': Subdivide integration domain in 2^floor(accuracy) points
%    - 'fastapproximation': Calculate the transmittance using an
%                              analytical approxmiation that evaluates faster. By default false.
%                              This approximation is very good for narrowband filters
%     - 'pixel' a pixel (see pixel2D) to change size of pixel relative to
%               filter size. By default the pixel will have the same size
%               as the filter. The pixel can not be outside of the filter
%   Outputs
%    - T (Wx1):  Ray-model estimation of the transmittance of a tiny Fabry-Pérot filter
%
%  Copyright Thomas Goossens
%  http://github.com/tgoossens/tinythinfilm
%

% % No inputparser because this makes things superslow
% %% Polarization recursion
% % If upolarized, recursively do the separate polarizations and average out
% % the transmittances
% if(or(polarization=='unpolarized',polarization=='unpolarised'))
%     [T_s] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidths,cwl,wavelengths,angledeg,'s',accuracy,pixelrange,flag_fastapproximation);
%     [T_p] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidths,cwl,wavelengths,angledeg,'p',accuracy,pixelrange,flag_fastapproximation);
%     T =  0.5*(T_s+T_p);
%     return;
% end

%% Variable argument saccuracy,flag_fastapproximation


%% Full thin film stack
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

% Transmission coefficient between incident medium and substrate (in absence of
% monolayer)
tsub =  2*eta0./(eta0+eta2);



% Pixel range
% The left of the filter by convention is at x=0; We need to do a coodinate
% transform. X=-filterwidth/2 becomes the origin
pixelwidth = pixelrange(2)-pixelrange(1);
pixel_start=pixelrange(1)-(-filterwidths(2)/2);
pixel_end=pixel_start+pixelwidth;


% Sampling of spatial pixel axis
x=linspace(pixel_start,pixel_end,2^floor(accuracy));
dx = abs(x(2)-x(1));



%% Right FILTER
% Calculation of number of reflections (anonymous fucntion)
% Filter dimensions
height_right=cwls(2)/(2*neff);
th_n = @(th) asind(sind(th)/neff);
num_right=@(x,th)x./(height_right*tand(th_n(th)));
N_right = @(x,th) floor(num_right(x,th)/2+1/2);
delta_right =(2*pi*neff*height_right*costh_n(2)./wavelengths') ; %



%% LEFT FILTER
% Calculation of number of reflections (anonymous fucntion)
% Filter dimensions
height_left=cwls(1)/(2*neff);
th_n = @(th) asind(sind(th)/neff);
num_left=@(x,th)x./(height_left*tand(th_n(th)));
N_left = @(x,th) floor(num_left(x,th)/2+1/2);

delta_left =(2*pi*neff*height_left*costh_n(2)./wavelengths') ; %



%% Reflectances
Rl = R(1);
Rr = R(2);






%% Refracted angle
theta_refracted =th_n(angledeg);

I_transmitted = zeros(numel(wavelengths),1);


%%
tl=sqrt(1-R(1));
tr=sqrt(1-R(2));
% Analytical approximation
% Numerical summation of the different contributions (with different
% number of interfering rays


for ix=1:numel(x)
    
    
    % Contributions of the main filter
    displacement_right = height_right*tand(theta_refracted);
    xi=x-displacement_right;
    gi = height_right*neff/cosd(theta_refracted);
    hi=height_right; %current filter
    
    E=0;
    n=0;
    while(xi>=filterwidths(1))
        di=(height_right-hi +xi*sind(angledeg))*cosd(angledeg);
        E = E+ t01*t12* (R(2)).^(n) .*exp(1i*(di+gi));
        
        xi=x_i-2*displacement_right; % (2 reflections)
        n=n+1;
    end
    
  

end
I_in = filterwidths(2); % We berekenen nog steeds de transmittanntie per pixel, dus alshet hoger ligt door cross-talk is dat ok en consistent.

% TODO aparte R voor elke filter
T =  (real(eta2)/real(eta0)).*I_transmitted/I_in;







end
