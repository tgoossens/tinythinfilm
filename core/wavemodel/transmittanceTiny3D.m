
function [T,Phi_t,Phi_in] = transmittanceTiny3D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixel)
% function [T,Phi_t,Phi_in] = transmittanceTiny3D(filter,incident_wavepacket,wavelengths,polarization,accuracy,pixelkernel)
%  transmittanceTiny3D  Simulate tiny filter transmittance
%    
%   Inputs
%    filter - Struct containing the tiny filter design (See also TINYFILTER)
%    incident_wavepacket -  A function @(nu,wavelength) of spatial frequency and wavelength that describes the angular spectrum of the field that enters the filter
%                           see directory wavepackets/
%    wavelengths (Wx1) - Wavelengths (same units as filter.width of filter)
%    polarization - ('s' or 'p' or 'unpolarized')    
%    accuracy - 2^floor(accuracy) subdivision of the spatial frequency domain.
%    pixelkernel - Encodes the width of the pixel and which spatial frequencies are sampled 
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - Phi_T (Wx1):  Transmitted flux [W]
%    - Phi_in (Wx1):  Incident flux [W]
%    
%  In the 3D case the spatial frequency nu is the amplitude on spatial
%  frequency = nu = sqrt(nu_x^2+nu_y^2)
%  
%  See also transmittanceTiny2D, tinyfilter, tinyfilter_equivalentmonolayer
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens

if(or(polarization=='unpolarized',polarization=='unpolarised'))
    [T_s,Phi_t_s,Phi_in_s] = transmittanceTiny3D(filter,incident_wavepacket,wavelengths,'s',accuracy,pixel);
    [T_p,Phi_t_p,Phi_in_p] =  transmittanceTiny3D(filter,incident_wavepacket,wavelengths,'p',accuracy,pixel);
    T =  0.5*(T_s+T_p);
    Phi_t =  0.5*(Phi_t_s+Phi_t_p);
    Phi_in =  0.5*(Phi_in_s+Phi_in_p);
    return;
end



wl=reshape(wavelengths,[1 1 numel(wavelengths)]); 

% Spatial frequency integration domain
nu_x = linspace(-1/wl(1), 1/wl(1),2^floor(accuracy))';
nu_y=nu_x';
nu = sqrt(nu_x.^2+nu_y.^2);
%nu = linspace(-10/wl(1), 10/wl(1),2^floor(accuracy))';


%% Definitions and helper ufcntions

% Wavenumber
k = @(n) 2*pi./(wl)*n;

% Fourier transform of the pixel kernel (so we don't recompute it for each wavelength)

pixelkernel_x = pixel.kernel.x(nu_x);
pixelkernel_y = pixel.kernel.y(nu_y);

%
%%  Calculate admittances

% Complex surface admittance of filter stack
% We will only use the transmission coefficient here
% Admittances of each layer


eta = admittance(filter.stack.refractiveindex(end),wl,nu,polarization);
eta_sub=eta(1);

eta = admittance(filter.stack.refractiveindex(1),wl,nu,polarization);
eta_in=eta(1);

%% Transmission coefficient
t = filter.transmission(wl,nu,polarization);

%% Simulate

for j=1:numel(wl)
    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wave
    Ain(:,:,j) = incident_wavepacket(filter.width,nu_x,nu_y,wl(j));
    
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

    function c = conv_pix(f)
        % If it only has zero elements, conv2fft will complain for some
        % reason (i.E. inplaceprod)
        numzero = sum(~(f(:)==0));
        if(numzero==0)
            c =  f;
            return
        end
        % The pixel kernel is separable so we can do 2 sequential 1D convolutions.
        c =conv2fft( conv2fft(f,pixelkernel_x,'same') ,pixelkernel_y,'same');
    end
end