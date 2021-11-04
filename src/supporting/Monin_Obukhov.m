function L = Monin_Obukhov(constants,meteo,H)

L           = -constants.rhoa*constants.cp*meteo.ustar.^3.* ...
            (meteo.Ta+273.15)./(constants.kappa*constants.g*H);           % [1]
%L(L<-1E3)   = -1E3;                                                     % [1] 
%L(L>1E2)    =  1E2;                                                     % [1]      
L(isnan(L)) = -1;
