function [T,Phi_t,Phi_in] = transmittanceTiny2DReflectiveBoundaries(filter,angledeg,wavelengths,polarization)
%  transmittanceTiny2DReflectiveBoundaries  Simulate tiny filter transmittance for perfectly reflecting (vertical) boundaries.  
%
%   [T, Phi_t, Phi_in] = transmittanceTiny2DReflectiveBoundaries(filter,angledeg,wavelengths,polarization,accuracy,pixelkernel);
%     
%   Inputs
%    - filter : Struct containing the tiny filter design (See also TINYFILTER)
%    - angledeg:  Incidence angle in degrees
%    - wavelengths (Wx1): Wavelengths (same units as filter.width of filter)
%    - polarization ('s' or 'p' or 'unpolarized')    
%    - accuracy: 2^floor(accuracy) subdivision of the spatial frequency domain.
%   Outputs
%    - T (Wx1):  Transmittance of the filter
%    - Phi_T (Wx1):  Transmitted flux [W]
%    - Phi_i (Wx1):  Incident flux [W]
%    
%    
%  See also TINYFILTER    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens

    if(or(polarization=='unpolarized',polarization=='unpolarised'))
        [T_s,Phi_t_s,Phi_in_s] =  transmittanceTiny2DReflectiveBoundaries(filter,angledeg,wavelengths,'s');
        [T_p,Phi_t_p,Phi_in_p] =   transmittanceTiny2DReflectiveBoundaries(filter,angledeg,wavelengths,'p');
        T =  0.5*(T_s+T_p);
        Phi_t =  0.5*(Phi_t_s+Phi_t_p);
        Phi_in =  0.5*(Phi_in_s+Phi_in_p);
        return;
    end
    
    
wl=reshape(wavelengths,[1 1 numel(wavelengths)]);
anglerad=deg2rad(angledeg); 
width=filter.width;


% Spatial frequency integration domain
% Exact evaluation at -1/wl(1) would result in a divison by zero for the p-polarized calculation

nu  = @(n) (n/(2*width));  % Q: why not n/(2*width)?


%%  Calculate admittances 

% Complex surface admittance of filter stack
% We will only use the transmission coefficient here

% Admittances of each layer

Phi_t = zeros(numel(wl),1);
Phi_in = zeros(numel(wl),1);


for j=1:numel(wl)
    

    nmax= floor(2*width/wl(j)); % Incidence angle 90 degrees.
    nrange=1:nmax; % No spatial frequencies beyond 90 degrees
        
    %%%%%%%%%% WAVE AMPLITUDES %%%%%%%%%%%%
    % Incident wwave
    k=2*pi./wl(j);
    A = @(n) 2*n.*pi.*(1-(-1).^n*exp(1i*k*width*sind(angledeg)))./(n.^2*pi^2-(k*width*sind(angledeg)).^2);    

    %[Y0,r,t] = surfaceadmittance(filter.stack.refractiveindex,filter.stack.thickness,wl(j),nu(nrange),polarization);
    t=filter.transmission(wl(j),nu(nrange),polarization);
    for n=1:numel(nrange)
        eta_sub = admittance(filter.stack.refractiveindex(end),wl(j),nu(nrange(n)),polarization);
        Phi_t(j) = Phi_t(j) + eta_sub(1)*abs(t(1,n))^2*abs(A(n))^2 * width;
            
    end

    Phi_t(j) = real(Phi_t(j))/4; % approximtaion sinus by exponentials (twice a factor of two-
    nu_angle=filter.stack.refractiveindex(1)*sin(anglerad)/wl(j);
    eta_in = admittance(filter.stack.refractiveindex,wl(j),nu_angle,polarization);
    Phi_in(j)=real(eta_in(1))/2 * filter.width;
    
    
end


    % Transmittance
    T=Phi_t./Phi_in;



end

