
function [T,Phi_t,Phi_in] = tinytransmittance3dfocus_rotational_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,coneangle_deg,cra_deg,alpha_deg,wavelengths,polarization,pixelkernel,accuracy);
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




if(or(polarization=='unpolarized',polarization=='unpolarised'))
    [T_s,Phi_t_s,Phi_in_s] = tinytransmittance3dfocus_rotational_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,coneangle_deg,cra_deg,alpha_deg,wavelengths,'s',pixelkernel,accuracy);
    [T_p,Phi_t_p,Phi_in_p] =  tinytransmittance3dfocus_rotational_mono(central_wavelength,normalized_fwhm,effective_index,width,nsub,coneangle_deg,cra_deg,alpha_deg,wavelengths,'p',pixelkernel,accuracy);
    T =  0.5*(T_s+T_p);
    Phi_t =  0.5*(Phi_t_s+Phi_t_p);
    Phi_in =  0.5*(Phi_t_s+Phi_t_p);
    return;
end


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
pixelkernel2D=pixelkernel_x.*pixelkernel_y;

%conv_pix=@(f) conv2fft(conv2fft(f,pixkernel_evaluated_x,'same'),pixkernel_evaluated_y,'same');

%
%
%% Definitions alpha and k
k = 2*pi./(wl)*neff;

alpha = @(v) sqrt(k.^2-(2*pi*v).^2);

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


%% Conversion functions
d=@(alpha) 0.5*pi*alpha; %delta
fwhm2r=@(alpha)-sqrt(cos(2*d(alpha)).^2-4*cos(2*d(alpha))+3)-cos(2*d(alpha))+2    ;

R = fwhm2r(normalized_fwhm);

delta = alpha(nu).*h;


ts=2*eta_in./(eta_sub+eta_in);

t = ts.*(1-R).*  1./(1-R*exp(1i*2*delta));



circ = @(u,radius) double(abs(u)<radius);



zi=1; % it doesn't matter for the pupil function (it is normalized out)
% however for comparison with the equations in article i keep it in

di=zi*tand(cra_deg);
lensradius=zi*tand(coneangle_deg);
P0 = @(x,y) circ(sqrt(x.^2+y.^2),lensradius);

P = @(x,y) P0(x+di*cosd(alpha_deg),y+di*sind(alpha_deg));





filteraperture = width^2 *sinca(pi*width*nu_x).*sinca(pi*width*nu_y);


% Allocate variables (not faster)
% Ain = zeros(numel(nu_x),numel(nu_y),numel(wl));
% At = zeros(numel(nu_x),numel(nu_y),numel(wl));
% Phi_in=zeros(numel(wl),1);
% Phi_t=zeros(numel(wl),1);

for j=1:numel(wl)
    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wwave
    width_y=width;
    

    pupilfunction= (P(-wl(j).*zi.*nu_x,-wl(j).*zi.*nu_y));
    %pupilfunction= (P(-wl.*zi.*nu_x,-wl.*zi.*nu_y)); % full vectorization
    %Ain(:,:,j) = conv2fft(pupilfunction,filteraperture,'same');
    Ain(:,:,j) = conv2fft(pupilfunction,filteraperture,'same');  
    
    % Useful integration domain;. This conditions corresponds to ignore incidence angles larger than 90 degres.
    domain = abs(nu).*wl(j) <=1;
    %domain = abs(nu).*wl <=1; 
    
    % Transmitted wave
    At(:,:,j)=domain.*t(:,:,j).*Ain(:,:,j);
    %At(:,:,:)=domain.*t(:,:,:).*Ain(:,:,:); 
    
    %%%%%%%%%% FLUXES  %%%%%%%%%%%%
    %Incident flux (explicit result)
    % nu_angle=sin(anglerad)/wl(j);
    % eta_in = admittance(1,wl(j),nu_angle,polarization);
    % Phi_in(j)=real(eta_in(1))/2 * width*width_y;
    
        % Transmitted flux
    temp=  0.5*real(eta_sub(:,:,j).*At(:,:,j).*conv_pix(conj(At(:,:,j))));
    temp= temp*abs(nu_x(2)-nu_y(1))*abs(nu_y(2)-nu_y(1)); % discretization convolution integral
    Phi_t(j)=trapz(nu_y,trapz(nu_x,temp,1),2);
      
    
    % Incident flux
    
    temp=  0.5*real(eta_in(:,:,j).*Ain(:,:,j).*conv_pix(conj(Ain(:,:,j))));
    temp= temp*abs(nu_x(2)-nu_y(1))*abs(nu_y(2)-nu_y(1)); % discretization convolution integral
    Phi_in(j)=trapz(nu_y,trapz(nu_x,temp,1),2);
    
%     % Transmitted flux
%     temp=  0.5*real(eta_sub(:,:,:).*At(:,:,:).*conv_pix(conj(At(:,:,:))));
%     temp= temp*abs(nu_x(2)-nu_y(1))*abs(nu_y(2)-nu_y(1)); % discretization convolution integral
%     Phi_t(:)=trapz(nu_y,trapz(nu_x,temp,1),2);
%       
%     
%     % Incident flux
%     
%     temp=  0.5*real(eta_in(:,:,:).*Ain(:,:,:).*conv_pix(conj(Ain(:,:,:))));
%     temp= temp*abs(nu_x(2)-nu_y(1))*abs(nu_y(2)-nu_y(1)); % discretization convolution integral
%     Phi_in(:)=trapz(nu_y,trapz(nu_x,temp,1),2);
    
end


% Transmittance
T=Phi_t./Phi_in;





    function f = sinca(x)
        % Modified sinc function because matlab sinc function already includes the factor pi.
        % This makes notation consistent with definitions in the publications.
        f=sinc(x/pi);
    end

    function c = conv_pix(f)
        
        %%
        % *conv2fft(f,pixelkernel2D,'same')* ;
        %conv2fft(conv2fft(f,pixelkernel_x,'same') ,pixelkernel_y);
        % efficient for large matrices
         %c =conv2fft(conv2fft(f,pixelkernel_x,'same') ,pixelkernel_y,'same');   ;
         c =conv2fft( conv2fft(f,pixelkernel_x,'same') ,pixelkernel_y,'same');
         % Separability makes it equally fast as using conv2
%        c = conv2(f,pixelkernel2D,'same');
        %c = conv2(c,pixelkernel_y,'same');  
        
    end
end