function  lens  = lensVignetted(exitpupil_distance,vignetcircle_radius,vignetcircle_sensitivity)
%anglesVignettedLens Summary of this function goes here
%   Detailed explanation goes here



% Return variable
lens.angulardistribution = @angularDistribution;
lens.name='Vignetted lens'


    function dA = angularDistribution(coneangle_deg,chiefray_deg,incidenceangle_deg)
        
        phideg=incidenceangle_deg;
        
        P = vignetcircle_radius;
        h = vignetcircle_sensitivity;
        x = exitpupil_distance;
        
        
        
        %% Circles
        eta = @(r,d,R) (real(acos((r.^2-R.^2+d.^2)./(2*d.*r)))); % finite aperture
        nu = @(r,d_r,P) (real(acos((r.^2-P.^2+d_r.^2)./(2*d_r.*r)))); % vignetting
        
        R=x*tand(coneangle_deg); % Radius of exit pupil
        
        % Displacements of the circles
        
        d   = x*tand(chiefray_deg);
        d_v = h*tand(chiefray_deg);
        d_r = d-d_v;
        
        %% Change of variable
        
        r=x*tand(phideg);  % radius of the blue circle in Fig. 4. & 5
        
        %% Calculate angular distribution
        gamma = NaN;
        if (h>=x)
            gamma = max(eta(r,d,R)-nu(r,d_r,P),0);
        elseif(h<x)
            gamma = min(eta(r,d,R),nu(r,d_r,P));
        end
        
        normalization=aperturearea(R,x,P,h,chiefray_deg);
        
        % Infinitesimal surface
        dA = 2*x^2*gamma*tand(phideg)./cosd(phideg).^2  *1./(normalization) ;
        
    end

    function A = aperturearea(R,x,P,h,cradeg)
        % APERTUREAREA Calculate the area of the exit pupil
        %
        % A = APERTUREAREA(R,x,P,h,cradeg). Chief ray angle in
        % degrees. Distances in mm.
        
        %
        % Asada, N., Amano, A., & Baba, M. (1996). Photometric calibration
        % of zoom lens systems. Proceedings - International Conference on
        % Pattern Recognition, 1,
        % 186–190. https://doi.org/10.1109/ICPR.1996.546016
        
        d   = x*tand(cradeg);
        d_v = h*tand(cradeg);
        
        % Helper functions
        alpha = real(acos((P^2-R^2+d_v^2)/(2*P*d_v)));
        beta = real(acos((R^2-P^2+d_v^2)/(2*R*d_v)));
        
        
        % area
        A = (P^2*(alpha-sin(alpha)*cos(alpha))+R^2*(beta-sin(beta)*cos(beta)));
        
    end

end


