clear all
clc
%% geometry-modified test version: no genetic regulation
tic
%% default coefficients (fixed)
n = 80; % number of cells for single layer embryo
R_y = 0.8; % raidus of yolk = radius of basal side of epithelium 
R_ap = 1.2; % radius of apical side of epithelium 
Ka = 10^4;   % cell area elasticity
Ky = 10^2; % yolk area elasticity
press = 10;
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
adhesion =10; % adhesion 
adsz = 15; % adhesion size
% predefined regions of tissue
tissue = 1:n;
vm = 0.5*tissue(1) + 0.25*tissue(end); % predefined ventral midline
center = [floor(vm) ceil(vm)]; % cell index at the midline
distance = tissue - vm; 
distance(distance>0.5*n) = distance(distance>0.5*n) - n;
mu = 0.02;

alpha0 = 9;
beta0 = 5.5;   
gamma0 = 1.0;  
baseline = 1;

alpha = alpha0*exp(-mu*distance.^2)+baseline;
beta = beta0 * ones(1,n);
gamma = gamma0 * ones(1,n);


[ap_l,bas_l,lat_l]=cal_length(pos);
ap_l_0 =mean(ap_l);


% ***********  iteration parameters ************** %
idx_frame = 1;
num_frame = 400;
N =  10^5;
ita = 1; % drag factor
dt = 1.5*10^(-4);
E = zeros(1,N);
for i = 1:N
%% energy minimization 
    
    [Egrad,Etot] =  cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,r0);
    net_force = -Egrad;
    E(i) = Etot;
    disp_f = dt/ita*net_force;
     pos = pos + disp_f;
    [ap_l,bas_l,lat_l]= cal_length(pos);




%% movie of invagination
if mod(i,num_frame)==1
%     H = axes;
%     drawfig(H,pos);
%     hold(H,'on')
%     plot(x_vm,y_vm,'--')
%     hold(H,'off')
%     p = get(gca,'position');
%     h = axes('parent',gcf,'position',[p(1)+.45 p(2)+.6 0.3*p(3) 0.2*p(4)]);
%     axis([0 80 0 12])
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
    
    
    Record.pos(idx_frame,:,:) = pos;
    Record.Aforce(idx_frame,:) = alpha;
    Record.Bforce(idx_frame,:) = beta;
    Record.Lforce(idx_frame,:) = gamma;
 
    idx_frame = idx_frame + 1;
end


    

 if mod(i,10000) == 0
     clc
     fprintf('progress: %d%% \n', i/N*100)
     drawfig(gca,pos)
     pause(0.01)
     toc
 end

end
% 
% H = axes;
% drawfig(H,pos);
% hold(H,'on')
% plot(x_vm,y_vm,'--')
% hold(H,'off')
% p = get(gca,'position');
% h = axes('parent',gcf,'position',[p(1)+.45 p(2)+.6 0.3*p(3) 0.2*p(4)]);
% axis([0 60 0 15])
% hold(h,'on')
% plot(h,alpha,'r')
% plot(h,beta,'b')
% plot(h,gamma,'g')
% pp = get(gca,'position');
% lgd = legend('apical','basal','lateral','location',[pp(1)+.14 pp(2)+.12 0.2*pp(3) 0.2*pp(4)]);
% legend('boxoff')
% title(lgd,'tensions')
% xlabel('cell index')
% hold(h,'off')



%movie(gcf,F,1,10,[0 0 0 0])




