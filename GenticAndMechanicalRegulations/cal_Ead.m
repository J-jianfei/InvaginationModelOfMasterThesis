function E_ad = cal_Ead(pos,invag_center,ref_l,ad_coef,ad_size)
%calculate cell-cell adhesion energy
n = length(pos)/2; % number of cells
inv_sz = 5; % check invagination curvature
ap_pos = pos(:,n+1:end);
global mutant 
if mutant == 2
    bas_pos = pos(:,1:n);
    if ~isempty(invag_center)
        m = size(invag_center,1);
        E = zeros(1,m);
        for i = 1:m
            
            p0 = invag_center(i,1) - inv_sz;
            p1 = invag_center(i,end);
            p2 = invag_center(i,end) + inv_sz;
            XData = bas_pos(1,[p0 p1 p2]);
            YData = bas_pos(2,[p0 p1 p2]);
            ta = sqrt((XData(2) - XData(1))^2 + (YData(2) - YData(1))^2);
            tb = sqrt((XData(3) - XData(2))^2 + (YData(3) - YData(2))^2);
            
            M(1,1) = 1;
            M(1,2) = -ta;
            M(1,3) = ta^2;
            M(2,1) = 1;
            M(2,2) = 0;
            M(2,3) = 0;
            M(3,1) = 1;
            M(3,2) = tb;
            M(3,3) = tb^2;
            
            a = M\XData';
            b = M\YData';
            
            kappa = abs(2*(a(3)*b(2) - b(3)*a(2))/(a(2)^2 + b(2)^2)^1.5);
            
            if kappa >= 1.25
                region_left = invag_center(i,1)-1 : -1:invag_center(i,1) - ad_size;
                region_right = invag_center(i,end)+2:invag_center(i,end)+ad_size+1;
                region_left(region_left<=0) = region_left(region_left<=0) + n;
                region_right(region_right>n) = region_right(region_right>n) - n;
                d = sqrt(sum((ap_pos(:,region_left) - ap_pos(:,region_right)).^2,1));
                E(i) = sum(ad_coef*d(d<ref_l));
            else
                E(i) = 0;
            end
            
            
        end
        E_ad = sum(E);
    else
        E_ad = 0;
    end
else
    if ~isempty(invag_center)
        for i = 1:size(invag_center,1)
            region_left = invag_center(i,1)-1 : -1:invag_center(i,1) - ad_size;
            region_right = invag_center(i,end)+2:invag_center(i,end)+ad_size+1;
            region_left(region_left<=0) = region_left(region_left<=0) + n;
            region_right(region_right>n) = region_right(region_right>n) - n;
            d = sqrt(sum((ap_pos(:,region_left) - ap_pos(:,region_right)).^2,1));
            E(i) = sum(ad_coef*d(d<ref_l));
        end
        E_ad = sum(E);
    else
        E_ad = 0;
    end
    
end


end

