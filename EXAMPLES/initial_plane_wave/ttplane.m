function [tP, tPs, tPpp, tPps, tPss] = ttplane(x0, xarray, p, H, Vp2, Vs2, Vp1, Vs1)
% this function computes the travel time of converted phases at
% Moho in a layer (Vp1, Vs1) over half space (Vp2, Vs2) model.
% [tP, tPs, tPpp, tPps, tPss] = ttplane(x0, xarray, p, H, Vp2, Vs2, Vp1, Vs1)
% the incoming plane wave is defined as f((x-x_0)*p+z*eta-t), where
% z is vertical pointing downwards and x is the horizontal direction.
% x can be a vector; H is in km, V is in km/s and p in s/km
% 
nx=length(xarray);
thp1=asin(p*Vp1); thp2=asin(p*Vp2);
ths1=asin(p*Vs1); ths2=asin(p*Vs2);
etap1=cos(thp1)/Vp1; etap2=cos(thp2)/Vp2;
etas1=cos(ths1)/Vs1; etas2=cos(ths2)/Vs2;
% initialization
tP=zeros(size(xarray)); tPs=tP; tPpp=tP; tPps=tP; tPss=tP;
for i = 1: nx
  x=xarray(i);
  if (x0+H*cot(thp2) > x) 
    error('Can not handle incidence from the right')
  end
  tt=p*(x-x0-H*cot(thp2));
  tP(i)=tt+etap1*H;
  tPs(i)=tt+etas1*H;
  tPpp(i)=tt+etap1*3*H;
  tPss(i)=tt+etas1*2*H+etap1*H;
  tPps(i)=tt+etap1*2*H+etas1*H;
% $$$   tP(i)=p*(x-xs)+etap1*H+etap2*zs;
% $$$   %Ps
% $$$   xx=xx1-H*tan(ths1);
% $$$   zs=xx*cos(thp2)*sin(thp2);
% $$$   xs=x2+xx*cos(thp2)^2;
% $$$   tPs(i)=p*(x-xs)+etas1*H+etap2*zs;
% $$$   %Ppp
% $$$   xx=xx1-2*H*tan(thp1);
% $$$   zs=xx*cos(thp2)*sin(thp2);
% $$$   xs=x2+xx*cos(thp2)^2;
% $$$   tPpp(i)=p*(x-xs)+etap1*2*H+etap2*zs;
% $$$   %Pps+Psp
% $$$   xx=xx1-H*tan(thp1)-H*tan(ths1);
% $$$   zs=xx*cos(thp2)*sin(thp2);
% $$$   xs=x2+xx*cos(thp2)^2;
% $$$   tPps(i)=p*(x-xs)+etap1*H+etas1*H+etap2*zs;
% $$$   %Pss
% $$$   xx=xx1-2*H*tan(ths1);
% $$$   zs=xx*cos(thp2)*sin(thp2);
% $$$   xs=x2+xx*cos(thp2)^2;
% $$$   tPss(i)=p*(x-xs)+etas1*2*H+etap2*zs;
end
  
  