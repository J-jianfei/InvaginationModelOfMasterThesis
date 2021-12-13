function [ap_l,bas_l,lat_l] = cal_length(pos)
% calculate apical/basal/lateral(connected to i-th vertex) length respectively, input:vertex positions
% note : lat_l is the sum of two lateral length of each cell
 num_c = length(pos)/2;
 ap_l = zeros(1,num_c);
 bas_l = zeros(1,num_c);

bas_v =pos(:,1:num_c);
bas_l(1:end-1) = sqrt(sum((bas_v(:,1:end-1) - bas_v(:,2:end)).^2,1));
bas_l(end) = sqrt(sum((bas_v(:,end)-bas_v(:,1)).^2,1));


ap_v = pos(:,num_c+1:end);
ap_l(1:end-1) = sqrt(sum((ap_v(:,1:end-1) - ap_v(:,2:end)).^2,1));
ap_l(end) = sqrt(sum((ap_v(:,end)-ap_v(:,1)).^2,1));

lat_l = sqrt(sum((ap_v - bas_v).^2,1));







end

