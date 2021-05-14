
function wavepacket =  wavepacket3DLensVignetted(coneangle_deg,cra_deg,azimuth_deg,vignetting_radius,vignetting_sensitivity,exitpupil_distance) 
% wavepacket3DLensVignetted
%  Wavepacket (angular spectrum) of light focused from a circular aperture (lens) onto a square filter of
%  width 'filterwidth'.
%
%  wavepacket = wavepacket3DLensVignetted(coneangle_deg,cra_deg,azimuth_deg,vignetting_radius,vignetting_sensitivity,exitpupil_distance) 
%     
%   Inputs
%     coneangle_deg - Half-cone angle of the focused light beam (related to
%                     f-number) in degrees
%     cra_deg - Chief ray angle in degrees
%     azimuth_deg - Off-axis direction
%     vignetting_radius - Radius of the vignetting circle
%     vignetting_sensitivity - such that the circle displacement d=vignetting_sensitivity*tand(cra_deg)
%     exitpupil_distance - Height of exit pupil above imaging plane
%
%   Outputs
%     wavepacket - An anonymous function @(filterwidth,nu_x,nu_y,wavelength) which
%                  represents the angular spectrum at specified wavelengths
%
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm
%    
    wavepacket = @field;
    function wavepacket_in = field(filterwidth,nu_x,nu_y,wavelength)   
        
        % Robustness feature, treat as rectangular if paired width
        % Treat as square filter if a single width is given
        if(numel(filterwidth)==2)
            filterwidth_x=filterwidth(1);
            filterwidth_y=filterwidth(2);
        else
            filterwidth_x=filterwidth;
            filterwidth_y=filterwidth_x;
        end


        % Height of pupil above imaging plane
        zi= exitpupil_distance; 
        
        % Distance of pixel from optical axis
        di=zi*tand(cra_deg);

        % Radius of the exit pupil given the half-cone angle and distance
        lensradius=zi*tand(coneangle_deg);
        
        
        % Circular aperture 
        P0 = @(x,y) circ(sqrt(x.^2+y.^2),lensradius);
        
        % Displaced vignetting circle 
        P0vignet = @(x,y) circ(sqrt(x.^2+y.^2),vignetting_radius);
        offset_vignettingcircle = -vignetting_sensitivity*tand(cra_deg); % dv
        Pvignet = @(x,y)P0vignet(x+offset_vignettingcircle*cosd(azimuth_deg),y+offset_vignettingcircle*sind(azimuth_deg));
        
        % The pupil function of the lens is the product of the two pupils,
        % this is equivalent to an INTERSECTION or logical AND operatio.
        P0lens = @(x,y) P0(x,y) .*Pvignet(x,y);
        pupilfunction = @(x,y) P0lens(x+di*cosd(azimuth_deg),y+di*sind(azimuth_deg));
    
        % Substitute the spatial frequencies (see derivation)
        pupilfunction_nu= (pupilfunction(-wavelength.*zi.*nu_x,-wavelength.*zi.*nu_y));
        
        % If the pupil function contains only zero elements, don't continue
        % because the convolution library will complain
        numzero = sum(~(pupilfunction_nu(:)==0));
        if(numzero==0)
            wavepacket_in =  pupilfunction_nu;
            return
        end
        
        % Apply the approximation of square aperture diffraction integral
        % This is modeled as a convolution operation (see Derivation)
        filteraperture = filterwidth_x.*filterwidth_y.*sinc_nopi(pi*filterwidth_x*nu_x).*sinc_nopi(pi*filterwidth_y*nu_y);
                
        wavepacket_in = conv2fft(pupilfunction_nu,filteraperture,'same');
    end


    
end