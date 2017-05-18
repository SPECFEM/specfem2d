function [ in ] = Inside_zone(x,z,zone, xout)
%
% [ in ] = Inside_zone(x,z,zone, xout)
%
%   test if (x,z) is inside zone
%   note (xout,zout) must be outside zone and zout is equal to z
%

in=0;
n=size(zone,1);
xtol=1e-3;
for i=1:1:n-1
    P0=zone(i,:);
    P1=zone(i+1,:);
    x0=P0(1);y0=P0(2);
    x1=P1(1);y1=P1(2);

    if (abs(x0 -x1) < xtol); 
        t=x0;
        if((t-x)*(t-xout) < 0 && (P0(2) - z)*(P1(2) - z) < 0);
            in=in+1;
        end
    else
        if ((P0(2) - z)*(P1(2) - z) < 0) 
            a=(y0-y1)/(x0-x1);
            b=(y1*x0-y0*x1) / (x0-x1);
            t = (z - b) / a;
            if ((t-x0)*(t-x1) < 0 && (t-x)*(t-xout) < 0);
                in=in+1;
            end
        end
    end

end



