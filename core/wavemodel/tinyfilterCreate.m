function filter = tinyfilterCreate(n_incident,n_stack,n_substrate,layerthickness,width)
%  TINYFILTER  Create a struct containing filter design parametrs
%   [filter] = TINYFILTER(n0,nstack,nsub,layerthickness,width)
%
%   Inputs
%    - n0: Refractive index of incident medium
%    - nsub: Refractive index of substrate
%    - n (Nx1) : Refractive indices of thin-film layers going from top layer to the layer on top of the subtrate.
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
%     filter.stack.refractiveindex=[1 1.5 2.4 1.5 1]
%     filter.stack.thickness=[NaN 10 20 10 NaN]  (NaN indicates incident and substrate medium (infinite thickness).
%     filter.width=5; %micron
%     filter.transmission = @(wavelength,spatialfrequency,polarization)
%     (...) where 'polarization'={'s','p','unpolarized'}
%    
%    
%    Copyright Thomas Goossens

    %% Checks
    assert(numel(layerthickness) == numel(n_stack));
    
    % Create struct
    filter.stack.refractiveindex=[n_incident reshape(n_stack,[1 numel(n_stack)]) n_substrate ];
    filter.stack.thickness=[NaN reshape(layerthickness,[1 numel(layerthickness)]) NaN ];
    filter.width=width;
    filter.transmission = @transmission;
    
    function t = transmission(wavelength,nu,polarization)
        [Y0,r,t] = surfaceadmittance(filter.stack.refractiveindex,filter.stack.thickness,wavelength,nu,polarization);
    end

end
