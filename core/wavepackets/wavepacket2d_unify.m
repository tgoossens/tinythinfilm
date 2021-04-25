 
function wavepacket =  wavepacket2d_collimated(incidence_angle_deg,incident_refractiveindex,filterwidth)
    
    wavepacket = @field;
    function wavepacket_in = field(nu_x,nu_y,wavelength)    
        nu=nu_x;
        wavepacket_in=filterwidth*sinca(pi*filterwidth*(nu-incident_refractiveindex*sind(incidence_angle_deg)./wavelength));
    end
    
    function f = sinca(x)
    % Modified sinc function because matlab sinc function already includes the factor pi.
    % This makes notation consistent with definitions in the publications.
        f=sinc(x/pi);
    end
    
end