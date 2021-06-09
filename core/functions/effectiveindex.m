function [neff] = effectiveindex(n_low,n_high,n_cavity,order)
% Effective refractive index for a given dielectric fabry perot filter.
%
%  INPUTS:
%     - n_low : Low refractive index 
%     - n_high: High refractive index 
%     - n_cavity: Cavity refractive index. This should be equal to either
%     n_low or n_high. The formula is different for both cases.
%     - order: Harmonic order (Main peak is order=1);
%
%   OUTPUTS:
%      - neff: Effective refractive index
% Origin of formulas: Macleod, H. A. (2017). Thin-Film Optical Filters,
% Fifth Edition. Thin-Film Optical Filters, Fifth Edition. CRC press. 
% Page 279
%
%
% Copyright Thomas Goossens

assert(order>=1);
m=order;


if(n_cavity == n_low)
    neff = n_low *sqrt((m-(m-1)*n_low/n_high)/(m-m*n_low/n_high+(n_low/n_high)^2));
elseif(n_cavity == n_high)
    neff = n_high * sqrt((m-(m-1)*n_low/n_high)./(m-1-(m-1)*n_low/n_high+n_high/n_low));
else
   warning('n_cavity should be either equal to n_low or n_high'); 
end

end

