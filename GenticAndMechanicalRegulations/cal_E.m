function E = cal_E(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0)
% calculate total energy of the epithelium: three line tension terms
% invovled regarding to Ana Hocevar Brezavscek's paper

n = length(pos)/2;
mat = [lat_l(2:length(lat_l)) lat_l(1)];
cell_area = cal_cell_area(pos);
yolk_area = polyarea(pos(1,1:n),pos(2,1:n));
E_cell = alpha.*ap_l + beta.*bas_l + 0.5*gamma.*(lat_l+mat)+0.5*Ka*(cell_area - ca0).^2;
E_yolk = 0.5*Ky*(yolk_area - ya0)^2;
E = sum(E_cell)+E_yolk;






end

