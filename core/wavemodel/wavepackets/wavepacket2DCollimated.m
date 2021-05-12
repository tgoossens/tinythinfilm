 
function wavepacket =  wavepacket2DCollimated(incidence_angle_deg,incident_refractiveindex)
% wavepacket2DCollimated
%  Wavepacket (angular spectrum) of a plane wave entering the filter at
%  oblique incidence
%
%  wavepacket =   wavepacket2DCollimated(incidence_angle_deg,incident_refractiveindex,filterwidth)
%     
%   Inputs
%     incidence_angle_deg - Incidence angle of the plane wave with respect
%                           to the normal vector of the surface.
%     incident_refractiveindex - Refractive index of the incident medium
%    
%
%   Outputs
%     wavepacket - An anonymous function @(filterwidth,nu,wavelength) which
%                  represents the angular spectrum at specified wavelengths
%
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm
    wavepacket = @field;
    function wavepacket_in = field(filterwidth,nu,wavelength)    
        filterwidth=filterwidth(1); % Use first dimension of potential 3D filter
        wavepacket_in=filterwidth*sinca(pi*filterwidth*(nu-incident_refractiveindex*sind(incidence_angle_deg)./wavelength));
    end
    
    function f = sinca(x)
    % Modified sinc function because matlab sinc function already includes the factor pi.
    % This makes notation consistent with definitions in the publications.
        f=sinc(x/pi);
    end
    
end