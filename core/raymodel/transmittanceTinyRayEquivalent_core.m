function [T,Tcoh] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,polarization,accuracy,pixelrange,flag_fastapproximation)
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
%    - filterwidth: Width of the film
%    - cwl: Central wavelength of the filter (same units as width)
%    - wavelengths (Wx1): Wavelengths (same units as height)
%    - angledeg:  Incidence angle in degrees
%    - polarization ('s' or 'p')
%    - 'accuracy': Subdivide integration domain in 2^floor(accuracy) points
%     - 'pixelrange' Range of the pixel (centered around the axis of the
%     filter). See also pixel2D for information.
%    - 'fastapproximation': Calculate the transmittance using an
%                              analytical approxmiation that evaluates faster. By default false.
%                              This approximation is very good for narrowband filters
%   Outputs
%    - T (Wx1):  Ray-model estimation of the transmittance of a tiny Fabry-Pérot filter
%
%  Copyright Thomas Goossens
%  http://github.com/tgoossens/tinythinfilm
%

% No inputparser because this makes things superslow
%% Polarization recursion
% If upolarized, recursively do the separate polarizations and average out
% the transmittances
if(or(polarization=='unpolarized',polarization=='unpolarised'))
    [T_s] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,'s',accuracy,pixelrange,flag_fastapproximation);
    [T_p] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,'p',accuracy,pixelrange,flag_fastapproximation);
    T =  0.5*(T_s+T_p);
    return;
end

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


% Filter dimensions
height=cwl/(2*neff);
w=filterwidth;

% Pixel range
% The left of the filter by convention is at x=0; We need to do a coodinate
% transform. X=-filterwidth/2 becomes the origin
pixelwidth = pixelrange(2)-pixelrange(1);
pixel_start=pixelrange(1)-(-filterwidth/2);
pixel_end=pixel_start+pixelwidth;


% Sampling of spatial pixel axis
x=linspace(pixel_start,pixel_end,2^floor(accuracy));
dx = abs(x(2)-x(1));

% Calculation of number of reflections (anonymous fucntion)
th_n = @(th) asind(sind(th)/neff);
num=@(x,th)x./(height*tand(th_n(th)));
N = @(x,th) floor(num(x,th)/2+1/2);

T = zeros(numel(wavelengths),1);
tcoh=T;
% Phase thickness

delta =(2*pi*equivstack(2)*height*costh_n(2)./wavelengths') ; %

%% Refracted angle
theta_refracted =th_n(angledeg);

% Analytical approximation
if(flag_fastapproximation)
    % Because of continuum approximation, we need to subtract pi, see supplementary material of Paper XX)
    delta_ct =delta- pi;
    
    % Maximum number of interfering rays, Bounded to avoid numerical
    % problems when M=infinity at normal incidence
    M1 = floor(min(pixel_start*neff/(cwl*tand(theta_refracted)),1e7));
    M2 = floor(min(pixel_end*neff/(cwl*tand(theta_refracted)),1e7));
    
    % Analytical equation
    denominator=((M2-M1).*(4*delta_ct.^2+log(R).^2));
    part2=-(R.^(2*M1)-R.^(2*M2))./(2*log(R).*(M2-M1));
    part3= (2*R.^(M1) .*(2*delta_ct.*sin(2*M1*delta)+cos(2*M1*delta_ct).*log(R)))./denominator;
    part4= -(2*R.^(M2) .*(2*delta_ct.*sin(2*M2*delta_ct)+cos(2*M2*delta_ct).*log(R)))./denominator;
    T=(real(eta2)/real(eta0))*(conj(tsub)*tsub).*(1-R).^2 .*(1+part2+part3+part4) ./(1-2*R*cos(2*delta_ct)+R.^2);;
    Tcoh=T; % DELETE
    
else
    % Numerical summation of the different contributions (with different
    % number of interfering rays
    for ix=1:numel(x)
        n = min(N(x(ix),angledeg),1e7);  % Number of interfering rays is limited to avoid division by zero when n=Inf
        formula=(R^(2*n)-2*R^(n)*cos(2*n*delta)+1)./(R^(2)-2*R*cos(2*delta)+1);
        

        part1=(R^n*exp(1i*2*n*delta)-1)./(R*exp(1i*2*delta)-1);
        
        tcoh = tcoh+dx*part1;
        T = T+dx*formula;
        
    end
    
    Tcoh =  (real(eta2)/real(eta0)) *(conj(tsub)*tsub)*(1-R)^2 * 1./filterwidth^2  *conj(tcoh).*tcoh   ;
    T =      (real(eta2)/real(eta0))*(conj(tsub)*tsub)*(1-R)^2*T/filterwidth;
end






end
