function Eext = cal_Eext(press,pos,ref1,ref2)
% calculate energy penalty moving outside vitelline
% if elliptic outline, ref1: length of major axis, ref2: length of minor
% axis
   % press = press*ones(1,length(pos));
   if nargin == 3 % hydrostatic global compression
       r = sqrt(sum(pos.^2,1));
       d = r - ref1;
       d(d<0) = 0;
       Eext = sum(press*(exp(d/ref1) - 1));
   elseif nargin == 4
       r = sqrt(sum(pos.^2,1));
       sin_theta = pos(2,:)./r;
       cos_theta = pos(1,:)./r;
       rho = sqrt(1./(cos_theta.^2/ref1^2 + sin_theta.^2/ref2^2));
       d = r - rho;
       d(d<0) = 0;
       Eext = sum(press*(exp(d./rho)-1));
   end

end

