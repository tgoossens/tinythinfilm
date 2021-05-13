function [T] = transmittanceTinyRayFocused(lens,n0,neff,nsub,R,width,cwl,wavelengths,coneangle_deg,chiefray_deg,polarization,accuracy,flag_fastapproximation) 
% transmittanceTinyRayEquivalent
% Simulate tiny transmittance of a Fabry-Pérot using a ray based model
% illuminted by focused light characterized by the half cone angle and
% chief ray angle.
%
% [T] = transmittanceTinyRayFocused(n0,neff,nsub,R,width,cwl,wavelengths,coneangle_deg,chiefray_deg,polarization,accuracy,flag_fastapproximation) 
%     
%   Inputs
%    - lens: a struct containing at least the field
%    'lens.angulardistribution  = @(coneangle,chiefray,incidenceangle)
%    - n0: Refractive index incident medium
%    - neff: effective refractive index of the cavity
%    - nsubstrate:Refractive index substrate
%    - R: Product of reflection coefficients
%    - width: Width of the film
%    - cwl: Central wavelength of the filter (same units as width)
%    - wavelengths (Wx1): Wavelengths (same units as height)
%    - coneangle_deg:  Half-cone angle in degrees (related to f-number)
%    - chiefray_deg:   Chief ray angle in degrees
%    - polarization ('s' or 'p')    
%    - accuracy: Subdivide integration domain in 2^floor(accuracy) points
%    - flag_fastapproximation: Calculate the transmittance using an
%    analytical approxmiation that evaluates faster. By default false. But
%    it works really well, you should try it ;)
%   Outputs
%    - T (Wx1):  Ray-model estimation of the transmittance of a tiny Fabry-Pérot filter 
%    
%  Copyright Thomas Goossens  
%  http://github.com/tgoossens/tinythinfilm
%


if ~exist('flag_fastapproximation','var')
    % Default value:
    flag_fastapproximation=false;;
end

%% 1D approximation angular distribution
conedeg = coneangle_deg;
cradeg=chiefray_deg;


%gamma = @(phideg) real(acos( (tand(cradeg+eps).^2-tand(conedeg).^2 +tand(phideg+eps).^2)./(2*tand(cradeg+eps).*tand(phideg+eps))));

%dA = anglesIdealLens(conedeg,cradeg);
dA = lens.angulardistribution;


% Sample to maximal angle in the distribution (according to geometrical
% optics).
phi_samples = linspace(0,ceil(atand(tand(cradeg)+tand(conedeg))),100); % radians

%%

for p=1:numel(phi_samples)
    temp(:,p)=dA(conedeg,cradeg,phi_samples(p))*transmittanceTinyRayEquivalent(n0,neff,nsub,R,width,cwl,wavelengths,phi_samples(p),polarization,accuracy,flag_fastapproximation);
end
temp(isnan(temp))=0;

% The derivaiton of the formula was in radians so the final integration
% needs to be done in the radian domain.
T= trapz(deg2rad(phi_samples),temp,2);



end
