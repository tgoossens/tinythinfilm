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

tl=sqrt(1-R(1));
tr=sqrt(1-R(2));
% Analytical approximation
% Numerical summation of the different contributions (with different
% number of interfering rays


for ix=1:numel(x)
    
    % Contribution of right filter (this is the main filter)
    Nr = min(N_right(x(ix),angledeg),1e7);  % Number of interfering rays is limited to avoid division by zero when n=Inf
    
    E_right =tr^2*tsub  *(Rr^Nr*exp(1i*2*Nr*delta_right)-1)./(Rr*exp(1i*2*delta_right)-1);
    
    % Contribution of left filter
    %%% 1. Get position in left filter: in coordinate system of the left
    %%% filter
    displacement_right = height_right*tand(theta_refracted);
    xleft = filterwidths(1)-(x(ix)-displacement_right*(2+2*(Nr-1)));
    %%%% 2. Apply ray model for the left filter
    Nl = min(N_left(xleft,angledeg),1e7);  % Number of interfering rays is limited to avoid division by zero when n=Inf
    
    %%%% 3. Determine the complex amplitude of ray arriving at this point
    E_left = tl*  (Rl.^Nl*exp(1i*2*Nl*delta_left)-1)./(Rl*exp(1i*2*delta_left)-1);
    
    % TODO: check of ik Rr en Rlcorect tel, want 1 reflectie is Rl, de
    % andere Rr
    % Het werkt veel beter als ik nog een scalar vermenidvuldigign doe met
    % *Rr^(Nr)  , waarom?
    E_cross_transmitted =  E_left .*((Rr)^Nr).^0 .*(Rr*exp(1i*2*delta_right)).^Nr * (tr*tsub);
    
    
    %%%%% INTEFEROMETRIC SUM OF ALL CONTRIBUTIONS
    
    % Je moet nog in rekening nemen dat rays andere fase hebben
    dphi=exp(-1i*2*pi*(height_left-height_right)/cosd(angledeg)./wavelengths');
    dphi=dphi+exp(-1i*2*pi*filterwidths(1)*sind(angledeg)./wavelengths');
    % this changes a lot!
    dphi=1;
    E_transmitted = E_right+E_cross_transmitted.*dphi;
    

    
    if(coherent_flag)
        % Inteferemoetric sym
        I_transmitted = I_transmitted + dx*conj(E_transmitted).*E_transmitted;
    else
        % Incoherent sum
        I_transmitted = I_transmitted+ dx*(conj(E_right).*E_right+conj(E_cross_transmitted).*E_cross_transmitted);
    end
    
end
I_in = filterwidths(2); % We berekenen nog steeds de transmittanntie per pixel, dus alshet hoger ligt door cross-talk is dat ok en consistent.

% TODO aparte R voor elke filter
T =  (real(eta2)/real(eta0)).*I_transmitted/I_in;







end
