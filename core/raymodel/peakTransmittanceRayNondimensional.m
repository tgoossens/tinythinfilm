function Tpeak = peakTransmittanceRayNondimensional(MlogR)
%peakTransmittanceRayNondimensional 

Tpeak = 1+ (1-exp(MlogR)).*(3-exp(MlogR))./(exp(2*MlogR));

end

