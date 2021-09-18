function [T,Tcoh] = transmittanceTinyRayEquivalent(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,polarization,varargin)
% transmittanceTinyRayEquivalent
% Simulate tiny transmittance of a Fabry-Pérot using a ray based model
%
% [T] = transmittanceTinyRayEquivalent(n0,neff,nsub,R,width,cwl,wavelengths,angledeg,polarization,accuracy)
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


%% Variable argument saccuracy,flag_fastapproximation
variableinputs = ieParamFormat(varargin);
p = inputParser;
p.addParameter('accuracy', 6, @isnumeric);
p.addParameter('fastapproximation',true,@islogical);
p.addParameter('pixel',NaN); % Set nan to avoid overhead of creating a pixel everytime

p.parse(variableinputs{:});


accuracy= p.Results.accuracy;
flag_fastapproximation= p.Results.fastapproximation;
pixel= p.Results.pixel;

if(~isstruct(pixel)) 
    pixel =struct;
    pixelrange=[-0.5 0.5]*filterwidth;
else 
    pixelrange=pixel.range.x;
end

if(or(polarization=='unpolarized',polarization=='unpolarised'))
    % Go immediately to the efficient implementation to avoid overhead of
    % inputparser
    [T_s,Tcoh_s] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,'s',accuracy,pixelrange,flag_fastapproximation);
    [T_p,Tcoh_p] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,'p',accuracy,pixelrange,flag_fastapproximation);
    T =  0.5*(T_s+T_p);
    Tcoh = 0.5*(Tcoh_s+Tcoh_p);
else
    
    [T,Tcoh] = transmittanceTinyRayEquivalent_core(n0,neff,nsub,R,filterwidth,cwl,wavelengths,angledeg,polarization,accuracy,pixelrange,flag_fastapproximation)
end


end
