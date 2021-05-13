function  lens  = lensIdeal()
%ANGLEDISTRIBUTIONIDEALLENS Summary of this function goes here
%   Detailed explanation goes here

lens.angulardistribution = @angularDistribution;
lens.name="Ideal Lens with circular aperture";


    function dA = angularDistribution(coneangle_deg,chiefray_deg,incidenceangle_deg)
        phideg=incidenceangle_deg;
        gamma =  real(acos( (tand(chiefray_deg+eps).^2-tand(coneangle_deg).^2 +tand(phideg+eps).^2)./(2*tand(chiefray_deg+eps).*tand(phideg+eps))));
        
        normalization= pi*tand(coneangle_deg).^2;
        dA = 2*gamma*tand(phideg)./cosd(phideg).^2  *1./(normalization) ;
        
    end

end

