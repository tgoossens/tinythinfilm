
function [T,Phi_t,Phi_in] = tinytransmittance3dfocus_core(filter,coneangle_deg,cra_deg,azimuth_deg,wavelengths,polarization,pixelkernel,accuracy);
%  TINYTRANSMITTANCE3D_MONO  Simulate tiny filter transmittance
%   [T] = TINYTRANSMITTANCE3D_MONO(filter,angledeg,wavelengths,polarization,accuracy);
%
%   Inputs
%    - filter : Struct containing the tiny filter design (See also TINYFILTER)
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



    %% TODO 
% if(or(polarization=='unpolarized',polarization=='unpolarised'))
%     [T_s,Phi_t_s,Phi_in_s] = tinytransmittance3dfocus_rotational_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,coneangle_deg,cra_deg,alpha_deg,wavelengths,'s',pixelkernel,accuracy);
%     [T_p,Phi_t_p,Phi_in_p] =  tinytransmittance3dfocus_rotational_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,coneangle_deg,cra_deg,alpha_deg,wavelengths,'p',pixelkernel,accuracy);
%     T =  0.5*(T_s+T_p);
%     Phi_t =  0.5*(Phi_t_s+Phi_t_p);
%     Phi_in =  0.5*(Phi_t_s+Phi_t_p);
%     return;
% end


wl=reshape(wavelengths,[1 1 numel(wavelengths)]); 
conerad=deg2rad(coneangle_deg);
crarad=deg2rad(cra_deg);

neff=effective_index;
cwl=central_wavelength;

% Spatial frequency integration domain
nu_x = linspace(-1/wl(1), 1/wl(1),2^floor(accuracy))';
nu_y=nu_x';
nu = sqrt(nu_x.^2+nu_y.^2);
%nu = linspace(-10/wl(1), 10/wl(1),2^floor(accuracy))';



%% Definitions and helper ufcntions

% Wavenumber
k = @(n) 2*pi./(wl)*n;

% Fourier transform of the pixel kernel (so we don't recompute it for each wavelength)
fftpix=fft(pixelkernel(nu_x));
pixelkernel_x = pixelkernel(nu_x);
pixelkernel_y = pixelkernel(nu_y);

%
%% Definitions azimuth and k
k = 2*pi./(wl)*neff;

azimuth = @(v) sqrt(k.^2-(2*pi*v).^2);

% Half wave plate
h = cwl/(2*neff);


%%  Calculate admittances

% Complex surface admittance of filter stack
% We will only use the transmission coefficient here
% Admittances of each layer

eta = admittance(nsub,wl,nu,polarization);
eta_sub=eta(1);

eta = admittance(1,wl,nu,polarization);
eta_in=eta(1);


%% Incident field
wavepacket_in=wavepacket_3dfocus(coneangle_deg,cra_deg,azimuth_deg,width);


%% Simulate

for j=1:numel(wl)
    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wwave

    Ain(:,:,j) = wavepacket_in(nu_x,nu_y,wl(j));
    
    % Useful integration domain;. This conditions corresponds to ignore incidence angles larger than 90 degres.
    domain = abs(nu).*wl(j) <=1;
    
    % Transmitted wave
    At(:,:,j)=domain.*t(:,:,j).*Ain(:,:,j);

    
    %%%%%%%%%% FLUXES  %%%%%%%%%%%%
    
        % Transmitted flux
    temp=  0.5*real(eta_sub(:,:,j).*At(:,:,j).*conv_pix(conj(At(:,:,j))));
    temp= temp*abs(nu_x(2)-nu_y(1))*abs(nu_y(2)-nu_y(1)); % discretization convolution integral
    Phi_t(j)=trapz(nu_y,trapz(nu_x,temp,1),2);
      
    % Incident flux
    temp=  0.5*real(eta_in(:,:,j).*Ain(:,:,j).*conv_pix(conj(Ain(:,:,j))));
    temp= temp*abs(nu_x(2)-nu_y(1))*abs(nu_y(2)-nu_y(1)); % discretization convolution integral
    Phi_in(j)=trapz(nu_y,trapz(nu_x,temp,1),2);
    
end


% Transmittance
T=Phi_t./Phi_in;



    function f = sinca(x)
        % Modified sinc function because matlab sinc function already includes the factor pi.
        % This makes notation consistent with definitions in the publications.
        f=sinc(x/pi);
    end

    function c = conv_pix(f)
    % The pixel kernel is separable so we can do 2 sequential 1D convolutions.
        c =conv2fft( conv2fft(f,pixelkernel_x,'same') ,pixelkernel_y,'same');
    end
end