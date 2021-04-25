
function wavepacket =  wavepacket3DLens(coneangle_deg,cra_deg,azimuth_deg,filterwidth) 
% wavepacket3DLens
%  Wavepacket (angular spectrum) of light focused from a circular aperture (lens) onto a square filter of
%  width 'filterwidth'.
%
%  wavepacket =   wavepacket3DLens(coneangle_deg,cra_deg,azimuth_deg,filterwidth)
%     
%   Inputs
%     coneangle_deg - Half-cone angle of the focused light beam (related to
%                     f-number) in degrees
%     cra_deg - Chief ray angle in degrees
%     azimuth_deg - Off-axis direction
%     filterwidth - Width in both dimensions of the square filter
%
%   Outputs
%     wavepacket - An anonymous function @(nu_x,nu_y,wavelength) which
%                  represents the angular spectrum at specified wavelengths
%
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm
%    
    wavepacket = @field;
    function wavepacket_in = field(nu_x,nu_y,wavelength)    
        circ = @(u,radius) double(abs(u)<radius);

        zi=1; % it doesn't matter for the pupil function (it is normalized out)
              % however for comparison with the equations in article i keep it in

        di=zi*tand(cra_deg);
        lensradius=zi*tand(coneangle_deg);
        P0 = @(x,y) circ(sqrt(x.^2+y.^2),lensradius);
        P = @(x,y) P0(x+di*cosd(azimuth_deg),y+di*sind(azimuth_deg));
        filteraperture = filterwidth^2 *sinca(pi*filterwidth*nu_x).*sinca(pi*filterwidth*nu_y);
        pupilfunction= (P(-wavelength.*zi.*nu_x,-wavelength.*zi.*nu_y));
        wavepacket_in = conv2fft(pupilfunction,filteraperture,'same');  
    end
    
    function f = sinca(x)
    % Modified sinc function because matlab sinc function already includes the factor pi.
    % This makes notation consistent with definitions in the publications.
        f=sinc(x/pi);
    end
    
end