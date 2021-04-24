function filter = tinyfilter_equivalentmonolayer(central_wavelength,normalized_fwhm,effective_index,width,refractiveindex_incident,refractiveindex_substrate)
%  TINYFILTER  Create a struct containing filter design parametrs
%   [filter] = TINYFILTER(n0,nstack,nsub,layerthickness,width)
%
%   Inputs
%    - n0: Refractive index of incident medium
%    - nsub: Refractive index of substrate
%    - n (Nx1) : Refractive indices of thin-film layers instack
%    - h (Nx1) : Thicknesses of thin-film layers in stack
%    - width: Width of the filter (same unit as h)
%
%   Outputs
%    - filter A struct conform the other functions in the library
%
%   Example
%    filter=tinyfilter(1,1,[1.5 2.4 1], [10 20 10],5)
%
%    Produces a struct equivalent to:
%     filter.n=[1 1.5 2.4 1.5 1]
%     filter.h=[NaN 10 20 10 NaN]  (NaN indicates incident and substrate medium (infinite thickness).
%     filter.width=5; %micron
%
%
%    Copyright Thomas Goossens

% Half wave plate thickness to center at the desired wavelength
h = central_wavelength/(2*effective_index);


filter.stack.refractiveindex=[refractiveindex_incident effective_index refractiveindex_substrate];
filter.stack.thickness=[NaN h NaN ];
filter.width=width;


filter.transmission = @transmission;

    function t = transmission(wavelength,nu,polarization)
        % Conversion functions
        d=@(alpha) 0.5*pi*alpha; %delta
        fwhm2r=@(alpha)-sqrt(cos(2*d(alpha)).^2-4*cos(2*d(alpha))+3)-cos(2*d(alpha))+2    ;

        % Wavenumber
        k = 2*pi./(wavelength)*effective_index;
        
        % Propagation coefficient
        alpha = @(v) sqrt(k.^2-(2*pi*v).^2);
        
        
        % Phase thickness
        delta = alpha(nu).*h;
        
        % Reflectance of the mirrors for a given normalized FWHM
        R = fwhm2r(normalized_fwhm);
        
        % Admittances
        eta = admittance(refractiveindex_substrate,wavelength,nu,polarization);
        eta_sub=eta(1);
        
        eta = admittance(refractiveindex_incident,wavelength,nu,polarization);
        eta_in=eta(1);
        
        
        % Transmission coefficient
        ts=2*eta_in./(eta_sub+eta_in);
        t = ts.*(1-R).*  1./(1-R*exp(1i*2*delta));
        
    end

end
