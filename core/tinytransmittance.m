function [T,Phi_t,Phi_in] = tinytransmittance(filter,angledeg,wavelengths,polarization,accuracy);
%  TINYTRANSMITTANCE  Simulate tiny filter transmittance
%   [T] = TINYTRANSMITTANCE(filter,angledeg,wavelengths,polarization,accuracy);
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
    
    
wl=wavelengths;
anglerad=deg2rad(angledeg);
width=filter.width;


nu = linspace(-1.05/wl(1), 1.05/wl(1),2^floor(accuracy))';
vv=nu;


dv = abs(vv(2)-vv(1));

%% Definitions alpha and k
k = @(n) 2*pi./(wl)*n;

alpha = @(v,n) sqrt(k(n).^2-(2*pi*v).^2);


%%  Calculate admittances
[Y,r,t] = surfaceadmittance(filter.n,filter.h,wl,nu,polarization);
eta = admittance(filter.n,wl,nu,polarization);

clear Eout;
%matlab sinc function already includes the factor pi!!!
%sinca = @(x) sin(x)./(x+eps);

sinca=@(x)sinc(x/pi);

% Create tilted wave front and do this per wavelength!
for j=1:numel(wl)
    
    Ain(:,j) = width*sinca(pi*width*(vv-filter.n(1)*sin(anglerad)/wl(j))); ...
    %custom
    bandpass = abs(vv).*wl(j) <=1;
    
    At(:,j)=bandpass.*t(:,j).*Ain(:,j);
    Etout(:,j)=fftshift(fft(At(:,j)));
end


Erout=zeros(size(Etout));

xpix = @(v) width*sinca(pi*width*v);
xfull=@(v) (v)==0;

%pix =xpix(vv'-vv);
fftpix=fft(xpix(vv));


% Admittances of incident and substrate medium
eta_in=eta(1);
eta_sub=eta(numel(filter.n));


for j=1:numel(wl)
    %OPGELET transpose in matlab neemt ook conj: gebruik .'

    %Incident flux
    temp = real(eta_in(:,j).*Ain(:,j).*fftshift(ifft(fft(conj(Ain(:,j))).*fftpix)));
    Phi_in(j)=trapz(vv,temp);

    
    % Transmutted flux
    Phi_t(j)=real(trapz(vv,eta_sub(:,j).*At(:,j).*fftshift(ifft(fft(conj(At(:,j))).* ...
                                                      fftpix))));
    
end

T=Phi_t./Phi_in;


end