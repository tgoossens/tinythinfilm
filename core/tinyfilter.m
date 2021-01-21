function filter = tinyfilter(n0,nstack,nsub,layerthickness,width)
%  TINYFILTER  Create a struct containing filter design parametrs
%   [filter] = SURFACEADMITTANCE(n,h,wl,v,polarization)
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
    
    
    filter.n=[n0 reshape(nstack,[1 numel(nstack)]) nsub ];
    filter.h=[NaN reshape(layerthickness,[1 numel(layerthickness)]) NaN ];
    filter.width=width;
    
    
end
