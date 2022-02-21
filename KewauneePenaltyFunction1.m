function [temp,penalty]=KewauneePenaltyFunction1(gammaG,alphaO,dcutoff,ncutoff,shapeval)

temp=((10000-alphaO)./10000).^(shapeval)+((10000-gammaG)./10000).^(shapeval); 


penalty=(temp/dcutoff).^(ncutoff); 




