 
function wavepacket =  wavepacket3DCollimated(incidence_angle_deg,azimuth_angle_deg,incident_refractiveindex)
% wavepacket3DCollimated
%  Wavepacket (angular spectrum) of an incidence plane wave on a filter of
%  width 'filterwidth'.
%
% wavepacket =   wavepacket3DCollimated(coneangle_deg,cra_deg,azimuth_deg,filterwidth)
%     
%   Inputs
%     incidence_angle_deg - Incidence angle w.r.t. to normal (zenith)
%     azimuth_angle_deg - Azimuth angle with respect to x-axis
%     incident_refractiveindex - Refractive index of incident medium
%     filterwidth - With of the filter (in both x and y)
%
%   Outputs
%     wavepacket - An anonymous function @(nu_x,nu_y,wavelength) which
%                  represents the angular spectrum at specified wavelengths
%
%  Coordinate system
%     (x,y) parallel with filter surface
%       z   perpendicular to the filter surface
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm
%



wavepacket = @field;

% The zenith and azimuth angle correspond to a linear slope in each
% horizontal dimension
angle_x = incidence_angle_deg*cosd(azimuth_angle_deg);
angle_y = incidence_angle_deg*sind(azimuth_angle_deg);

    function wavepacket_in = field(filterwidth,nu_x,nu_y,wavelength)
        if(numel(filterwidth)==2)
            filterwidth_x=filterwidth(1);
            filterwidth_y=filterwidth(2);
        else
            filterwidth_x=filterwidth;
            filterwidth_y=filterwidth_x;
        end
        
        % We can make use of separability
        wavepacket_in_x=filterwidth_x*sinca(pi*filterwidth_x*(nu_x-incident_refractiveindex*sind(angle_x)./wavelength));
        wavepacket_in_y=filterwidth_y*sinca(pi*filterwidth_y*(nu_y-incident_refractiveindex*sind(angle_y)./wavelength));
        wavepacket_in = wavepacket_in_x .* wavepacket_in_y;
    end

    function f = sinca(x)
        % Modified sinc function because matlab sinc function already includes the factor pi.
        % This makes notation consistent with definitions in the publications.
        f=sinc(x/pi);
    end

end