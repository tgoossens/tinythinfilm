function  dA  = anglesVignettedLens(coneangle_deg,chiefray_deg)
%ANGLEDISTRIBUTIONIDEALLENS Summary of this function goes here
%   Detailed explanation goes here

gamma = @(phideg) real(acos( (tand(chiefray_deg+eps).^2-tand(coneangle_deg).^2 +tand(phideg+eps).^2)./(2*tand(chiefray_deg+eps).*tand(phideg+eps))));

normalization= pi*tand(conedeg).^2;
dA = @(phideg) 2*gamma(phideg)*tand(phideg)./cosd(phideg).^2  *1./(normalization) ;


end

