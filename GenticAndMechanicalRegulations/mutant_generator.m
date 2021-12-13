clear all
clc
% global R_y alpha beta gamma ya0 ca0 Ka Ky
n = 80; % number of cells for single layer embryo
R_y = 0.8; % raidus of yolk = radius of basal side of epithelium 
R_ap = 1.2; % radius of apical side of epithelium 
Ka = 10000;   % cell area elasticity
Ky = 100; % yolk area elasticity
idx = 1:n;
x_bas = R_y*cos((idx - 1)*2*pi/n );
y_bas = R_y*sin((idx - 1)*2*pi/n);
x_apc = R_ap*cos((idx - 1)*2*pi/n);
y_apc = R_ap*sin((idx - 1)*2*pi/n);
pos = [x_bas x_apc;y_bas y_apc];
yolk_area = polyarea(x_bas,y_bas);
ya0 = yolk_area; % prefered yolk areads
cell_area = cal_cell_area(pos);
ca0 = mean(cell_area); % prefered cell area
x_vm = [pos(1,n+1:end) pos(1,n+1)]; 
y_vm = [pos(2,n+1:end) pos(2,n+1)]; % get vitelline membrane position
r0 = mean(sqrt(x_vm(1:end-1).^2 + y_vm(1:end-1).^2)); % radius of membrane
press = 10;

adhesion = 10; % adhesion 
adsz = 15; % adhesion size
% predefine regions of tissue
tissue = 1:n;
% vm = 0.5*tissue(1) + 0.25*tissue(end); % predefined ventral midline
vm = 0.5*tissue(1) + 0.75*tissue(end); % for ectopic compression
center = [floor(vm) ceil(vm)]; % cell index at the midline
distance = tissue - vm; 
distance(distance>0.5*n) = distance(distance>0.5*n) - n;

[ap_l,bas_l,lat_l]=cal_length(pos);
ap_l_0 =mean(ap_l);

mu = 0.025;
alpha0 = 5;
beta0 = 5.0;
gamma0 = 1.0;
baseline = 1;

alpha = alpha0*exp(-mu*distance.^2)+baseline;
beta = beta0 * ones(1,n);
gamma = gamma0 * ones(1,n);

%%% uncomment: ectopic compression
% mid = 60.5;
% ecto_region = [36:n 1:15];
% ecto_region2 = 36:95;
% a_grad = 0.2;
% alpha = ones(1,n);
% alpha(ecto_region(ecto_region2<=mid)) = alpha(ecto_region(ecto_region2<=mid)) + a_grad*(ecto_region2(ecto_region2<=mid)-ecto_region(1));
% alpha(ecto_region(ecto_region2>floor(mid))) = alpha(ecto_region(ecto_region2>floor(mid))) - a_grad*(ecto_region2(ecto_region2>floor(mid))-ecto_region2(end));
% beta = ones(1,n); 
% gamma = ones(1,n);
%% uncomment: gloabal compression

% A0 = 5; B0 = 1; C0 = 1;
% alpha = A0*ones(1,n);
% beta = B0*ones(1,n);
% gamma = C0*ones(1,n);

%% elliptic compression
% alpha = 0.1*ones(1,n);
% beta = 2*ones(1,n);
% gamma = 0.1*ones(1,n);


idx_frame = 1;
num_frame = 150;
N =  0.5*10^4;
ita = 1; % drag factor
dt = 1.5*10^(-4);


for i = 1:N
    
    %  Etot0 = cal_Etot(ap_l,bas_l,lat_l,pos,center,0.5*ap_l_0,adhesion,adsz,r0,press);
   net_force = - cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,1*r0,1*r0);
    disp_f = dt/ita*net_force;
    pos = pos + disp_f;
       [ap_l,bas_l,lat_l]= cal_length(pos);
    
    
    
%     if mod(i,num_frame)==1
%     H = axes;
%     drawfig(H,pos);
% %     hold(H,'on')
% %     plot(x_vm,y_vm,'--')
% %     hold(H,'off')
%     p = get(gca,'position');
%     h = axes('parent',gcf,'position',[p(1)+.45 p(2)+.6 0.3*p(3) 0.2*p(4)]);
%     axis([0 60 0 12])
%     hold(h,'on')
%     plot(h,alpha,'r')
%     plot(h,beta,'b')
%     plot(h,gamma,'g')
%     pp = get(gca,'position');
%     lgd = legend('apical','basal','lateral','location',[pp(1)+.14 pp(2)+.12 0.2*pp(3) 0.2*pp(4)]);
%     legend('boxoff')
%     title(lgd,'tensions')
%     xlabel('cell index')
%     hold(h,'off')   
%     F(idx_frame) = getframe(gcf);       
%     idx_frame = idx_frame + 1;
%     end
 if mod(i,10000) == 0
     fprintf('progress: %d percents\n', i/N*100)
 end

end
    
    
H = axes;
drawfig(H,pos);
hold(H,'on')
plot(x_vm,y_vm,'--')
hold(H,'off')
p = get(gca,'position');
h = axes('parent',gcf,'position',[p(1)+.45 p(2)+.6 0.3*p(3) 0.2*p(4)]);
axis([0 80 0 12])
hold(h,'on')
plot(h,alpha,'r')
plot(h,beta,'b')
plot(h,gamma,'g')
pp = get(gca,'position');
lgd = legend('apical','basal','lateral','location',[pp(1)+.14 pp(2)+.12 0.2*pp(3) 0.2*pp(4)]);
legend('boxoff')
title(lgd,'tensions')
xlabel('cell index')
hold(h,'off')
%     
    
    
    
    
    
    
    
    
    
    
    
    