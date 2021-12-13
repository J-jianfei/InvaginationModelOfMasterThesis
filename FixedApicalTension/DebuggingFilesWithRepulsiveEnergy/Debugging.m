clear all
clc

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
beta0 = 1.5;
gamma0 = 0.5;
baseline = 1;

alpha = alpha0*exp(-mu*distance.^2)+baseline;
beta = beta0 * ones(1,n);
gamma = gamma0 * ones(1,n);


[ap_l,bas_l,lat_l]=cal_length(pos);
ap_l_0 =mean(ap_l);


% ***********  iteration parameters ************** %
idx_frame = 1;
num_frame = 300;
N =  10^5;
ita = 1; % drag factor
dt = 1.5*10^(-4);
E = zeros(1,N);

% load('bug1.mat','pos');
% load('Demo.mat')
load('bug3.mat')
pos = squeeze(Record.pos(232,:,:));
[ap_l,bas_l,lat_l]= cal_length(pos);
%%
for i = 1:400
    [Egrad,Etot] =  cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,r0);
    net_force = -Egrad;
    disp_f = dt/ita*net_force;
    pos = pos + disp_f;
    [ap_l,bas_l,lat_l]= cal_length(pos);
    
    drawfig(gca,pos)
    ff(i) = getframe(gcf);
%     if i == 315
%         pause
%     end
%     
end
% drawfig(gca,pos)