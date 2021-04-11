function [T,Phi_t,Phi_in] = tinytransmittance_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,polarization,accuracy)
%  TINYTRANSMITTANCE_MONO  Simulate tiny filter transmittance for an equivalent monolayer model
%   [T] = TINYTRANSMITTANCE_MONO(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,polarization,accuracy);
%    
%   Inputs
%    - central_wavelength: Central wavelength of the filter
%    - normalized_fwhm:  Lambda = FWHM / central_wavelength    
%    - effective_index: Effective refractive index 
%    - width: Width of the filter
%    - nsub: Substrate refractive index
%    - angledeg:  Incidence angle in degrees
%    - wavelengths (Wx1): Wavelengths (same units as filter.width of filter)
%    - polarization ('s' or 'p')    
%    - accuracy: 2^floor(accuracy) subdivision of the spatial frequency domain.
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - Phi_T (Wx1):  Transmitted flux [W]
%    - Phi_T (Wx1):  Incident flux [W]
%    
%    
%  See also TINYFILTER    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens
    

    if(or(polarization=='unpolarized',polarization=='unpolarised'))
        [T_s,Phi_t_s,Phi_in_s] = tinytransmittance_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,'s',accuracy);
        [T_p,Phi_t_p,Phi_in_p] = tinytransmittance_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,angledeg,wavelengths,'p',accuracy);
        T =  0.5*(T_s+T_p);
        Phi_t =  0.5*(Phi_t_s+Phi_t_p);
        Phi_in=  0.5*(Phi_in_s+Phi_in_p);
        return;
    end
    
    
    % Input variables
    wl=reshape(wavelengths,[1 1 numel(wavelengths)]);
    anglerad=deg2rad(angledeg); 
    neff=effective_index;
    cwl=central_wavelength;

    % Spatial frequency integration domain
    nu = linspace(-1/wl(1), 1/wl(1),2^floor(accuracy))';


    %% Definitions and helper ufcntions

    % Wavenumber
    k = @(n) 2*pi./(wl)*n; 

    % Fourier transform of the pixel kernel (so we don't recompute it for each wavelength)
    pixelkernel=pixel_fullwidth(width);
    pixelkernel_eval = pixelkernel(nu);
    conv_pix=@(f) conv(f,pixelkernel(nu),'same');


    %% Definitions alpha and k

    k = 2*pi./(wl)*neff;
    alpha = @(v) sqrt(k.^2-(2*pi*v).^2);

    % Half wave plate
    h = cwl/(2*neff);

    %%  Calculate admittances

    eta = admittance(nsub,wl,nu,polarization);
    eta_sub=eta(1);


    eta = admittance(1,wl,nu,polarization);
    eta_in=eta(1);



    %% Convert normalized FWHM to mirror reflectance

    d=@(alpha) 0.5*pi*alpha; %delta
    fwhm2r=@(alpha)-sqrt(cos(2*d(alpha)).^2-4*cos(2*d(alpha))+3)-cos(2*d(alpha))+2    ;
    R = fwhm2r(normalized_fwhm);



    delta = alpha(nu).*h;

    %% Convert to 2D matrices
    % The core functions put wavelength in the 3th dimension to faciliate extensions to 3D. However can reshape the matrix to have only rows and columns such that it can be processed by MATLAB's convolution functions

    eta_in=reshape(eta_in,[size(eta_in,1) size(eta_in,3)]);
    eta_sub=reshape(eta_sub,[size(eta_sub,1) size(eta_sub,3)]);
    delta = reshape(delta,[size(delta,1) size(delta,3)]);


    %% Transmission coefficient
    ts=2*eta_in./(eta_sub+eta_in); 
    t = ts.*(1-R).*  1./(1-R*exp(1i*2*delta));

    % reshape 2 normal 2D matrix, this is is because the code relies on the 
    % generalizeation to 3D
    wl2d=reshape(wl,[1 numel(wl)]);

    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wwave
    Ain = width*sinca(pi*width*(nu-sin(anglerad)./wl2d)); 
    
    % Useful integration domain;. This conditions corresponds to ignore incidence angles larger than 90 degres.
    domain = abs(nu).*wl2d <=1;
    
    % Transmitted wave
    At=domain.*t.*Ain;
    
    %%%%%%%%%% FLUXES  %%%%%%%%%%%%
    %Incident flux (explicit result)
    nu_angle=sin(anglerad)./wl;
    eta_in = admittance(1,wl,nu_angle,polarization);
    Phi_in=(real(eta_in(1))/2 * width);
    
    % Transmitted flux
    %    temp=  0.5*real(eta_sub(:,1,j).*At(:,1,j).*conv_pix(conj(At(:,1,j))));
    temp=  0.5*real(eta_sub.*At.*conv2fft(conj(At(:,:)),pixelkernel_eval,'same'));
    temp= temp*abs(nu(2)-nu(1)); % discretization convolution integral
    Phi_t=(trapz(nu,temp)); 

    % Transmittance
    T=Phi_t(:)./Phi_in(:);

    function f = sinca(x)
    % Modified sinc function because matlab sinc function already includes the factor pi.
    % This makes notation consistent with definitions in the publications.    
        f=sinc(x/pi); 
    end
end

