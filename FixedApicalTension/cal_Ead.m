function E_ad = cal_Ead(pos,invag_center,ref_l,ad_coef,ad_size)
%calculate cell-cell adhesion energy
n = length(pos)/2; % number of cells
ap_pos = pos(:,n+1:end);
% bas_pos = pos(:,1:n);
if ~isempty(invag_center)
    region_left = invag_center(1)-1 : -1:invag_center(1) - ad_size;
    region_right = invag_center(end)+2:invag_center(end)+ad_size+1;
    region_left(region_left<=0) = region_left(region_left<=0) + n;
    region_right(region_right>n) = region_right(region_right>n) - n;
    d = sqrt(sum((ap_pos(:,region_left) - ap_pos(:,region_right)).^2,1));
    E_ad = sum(ad_coef*d(d<ref_l));
  %  E_ad = sum(ad_coef*(1./(exp(d/ref_l)+1)));
elseif isempty(invag_center)
    E_ad = 0;
end

end

