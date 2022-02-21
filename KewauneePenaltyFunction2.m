function penalty=KewauneePenaltyFunction2(alphaO,dcutoff,ncutoff,shapeval)

temp=((10000-alphaO)./10000).^(shapeval);

penalty=(temp/dcutoff).^ncutoff;