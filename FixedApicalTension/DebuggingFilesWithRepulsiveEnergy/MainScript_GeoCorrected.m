 clear all
clc
%% combined-reponse version: high mechanical contribution
tic
%% default coefficients (fixed)
% global R_y alpha beta gamma ya0 ca0 Ka Ky
n = 60; % number of cells for single layer embryo
R_y = 1; % raidus of yolk = radius of basal side of epithelium 
R_ap = 1.2; % radius of apical side of epithelium 
Ka = 10000;   % cell area elasticity
Ky = 100; % yolk area elasticity
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
adhesion = 10; % adhesion 
adsz = 15; % adhesion size
% predefine regions of tissue
tissue = 1:n;
vm = 0.5*tissue(1) + 0.25*tissue(end); % predefined ventral midline
center = [floor(vm) ceil(vm)]; % cell index at the midline
distance = tissue - vm; 
distance(distance>0.5*n) = distance(distance>0.5*n) - n;

[ap_l,bas_l,lat_l]=cal_length(pos);
ap_l_0 =mean(ap_l);


% ***********  iteration parameters ************** %
idx_frame = 1;
num_frame = 400;
N =  10^5;
ita = 1; % drag factor
dt = 1.5*10^(-4);

%% parameters of genetic regulation

r1 = 0.75; % relocalization coeff from basal to apical  default:0.75

ut = 0.8; % twist breakdown rate;
um = 1; % myosin breakdown rate;
um_l = 1.2;
pt = 1; % twist production rate;
pm = 1; % constant myosin production rate;
%p0 = 1.2; % constant twist production rate modulated by a Gaussian to mimic the position dependence
ka = 1; % saturation of twist  % defaultL: 1
lm = 1.25; % basal myosin loss rate on mesoderm %  default:1.5
pm_a = 1; % apical myosin concentration
pm_e = 0.1; % basal myosin production rate
um_e = 0.01; % basal myosin breakdown rate 
mu = 0.02; % std of Gaussian term
kh =1;
pm_l = 2; % lateral constant myosin production 2
Th = 0.5; % reference value in Hill term
q = 4; % power in hill type mechanical response
p_max = 1.2; % maximal amplitude of Hill term
a = 0.7;
b = 3.75; %3.75


fb = 0.1; % basal myosin concentration to tension factor  0.2
fa = 1.2; % apical myosin concentration to tension factor
fg = 1; % lateral myosin concentration to tension factor

%trigger = zeros(n,1); % mechanical response as trigger 


% % twist concentration profile over all cells 

cm_apc = zeros(n,N);
cm_apc(:,1) = 0.0;
s = 0.001; % a constant to prevent from 0 denominator

% generate mutant with initial twist gradient and 
    ct = zeros(n,N);
    mutant = 0; 
    if mutant == 1
%         %     c0  = 1.5;
%         %     ct(:,1) = c0*exp(-mu*distance'.^2);
%         x_bas_init = 0.85*R_y*cos((idx - 1)*2*pi/n );
%         y_bas_init = 0.85*R_y*sin((idx - 1)*2*pi/n);
%         x_apc_init = 0.85*R_ap*cos((idx - 1)*2*pi/n);
%         y_apc_init = 0.85*R_ap*sin((idx - 1)*2*pi/n);
%         r_mut = mean(sqrt(x_apc_init.^2 + y_apc_init.^2));
%         press_mut = 100;
%         pos =[x_bas_init x_apc_init;y_bas_init y_apc_init];
%         [ap_l,bas_l,lat_l] = cal_length(pos);
%         applied_time = 1*N;
           pos = load('pos_ellipse.mat').pos;
           [ap_l,bas_l,lat_l] = cal_length(pos);
           x_vm = pos(1,n+1:end);
           y_vm = pos(2,n+1:end);
    elseif mutant == 2
        
% ********* keep the shape **********%
        pos = load('pos_ecto_compressed.mat').pos;
        [ap_l,bas_l,lat_l] = cal_length(pos);
% %         applied_time = 1;
    elseif mutant == 3
        % ********* keep the force ****** %
        pos = load('pos_ecto_compressed.mat').pos;
        [ap_l,bas_l,lat_l] = cal_length(pos);
        applied_time = N;  % number of iterations to keep the ecto compressive force
        mid = 45.5;
        ecto_region = [26:n 1:5];
        ecto_region2 = 26:65;
        a_grad = 0.4;
        alpha_ext = zeros(1,n);
        alpha_ext(ecto_region(ecto_region2<=mid)) = alpha_ext(ecto_region(ecto_region2<=mid)) + a_grad*(ecto_region2(ecto_region2<=mid)-ecto_region(1));
        alpha_ext(ecto_region(ecto_region2>floor(mid))) = alpha_ext(ecto_region(ecto_region2>floor(mid))) - a_grad*(ecto_region2(ecto_region2>floor(mid))-ecto_region2(end));
%         beta = ones(1,n);
%         gamma = ones(1,n);
         
    end

   

cm_bas = zeros(n,N);
cm_bas(:,1) = 30;
cm_lat =zeros(n,N);

alpha_min = 0.1;
beta_min = 0.1;
gamma_min = 0.1;


% r_inner = zeros(1,N); % inner radius
% r_min = zeros(1,N); % minimum radius during invagination
% r_16 = zeros(1,N); % radius of invagination center
% sp = 1;
% j = 1;
% time_interval = 100;

Record = struct('pos',[],'Aforce',[],'Bforce',[],'Lforce',[]);
for i = 1:N
%% translation concentration to line tension


%% energy minimization 

alpha = fa*cm_apc(:,i)';
beta = fb*cm_bas(:,i)';
gamma = fg*cm_lat(:,i)';
% gamma = 2*ones(1,n);
alpha(alpha<alpha_min) = alpha_min;
beta(beta<beta_min) = beta_min;
gamma(gamma<gamma_min) = gamma_min;


if mutant == 0    
    net_force = - cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,r0);
elseif mutant == 1 % elliptic global compression        
    net_force = - cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,0.9*r0,1.1*r0);
elseif mutant == 2       % local compression, keep the shape
    if i >= 0.3*N
        net_force = - cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,r0);
    else
        net_force = 0;
    end  
elseif mutant == 3  % local compression, keep the force
    if i <= applied_time
        alpha = alpha + alpha_ext;
    end
       net_force = - cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,r0); 
end



disp_f = dt/ita*net_force;
pos = pos+disp_f;
[ap_l,bas_l,lat_l]= cal_length(pos);

%% regulation network
% ********* try to add combined mechanical and gentic response ******* %
% L = (ap_l_0-ap_l')/ap_l_0;
% L(L<0) = 0;
% dct = dt*(pt*ct(:,i).^2./(1+ka*ct(:,i).^2) - ut*ct(:,i) + p_max*(a*exp(-mu*distance'.^2) + b*L).^q./(Th^q + (a*exp(-mu*distance'.^2) + b*L).^q));
L = abs(ap_l_0 - ap_l')/ap_l_0;
% L(L<=0.05) = 0;
prod = p_max*(a*exp(-mu*distance'.^2) + b*L).^q./(Th^q + (a*exp(-mu*distance'.^2) + b*L).^q);
dct = dt*(pt*ct(:,i).^2./(1+ka*ct(:,i).^2) - ut*ct(:,i) + prod);
dcm_apc = dt*(r1*ct(:,i).^2- um*cm_apc(:,i)+pm_a);
dcm_bas = dt * ( - lm*ct(:,i).^2 + pm_e*cm_bas(:,i) - um_e*cm_bas(:,i).^2);
 dcm_lat = dt * ( (lm - r1)*ct(:,i).^2./(kh*(cm_apc(:,i)+s)) - um*cm_lat(:,i) + pm_l); 
%    dcm_lat = dt * ( (lm - r1)*ct(:,i).^2 - um_l*cm_lat(:,i) + pm_l);
    cm_lat(:,i+1) = cm_lat(:,i) + dcm_lat;
ct(:,i+1) = ct(:,i) + dct;
ct(ct(:,i+1)<=0,i+1) = 0 ;
cm_apc(:,i+1) = cm_apc(:,i)+dcm_apc;
cm_bas(:,i+1) = cm_bas(:,i) + dcm_bas;

%% movie of invagination
if mod(i,num_frame)==1
            H = axes;
            drawfig(H,pos);
            hold(H,'on')
            plot(x_vm,y_vm,'r--')
 %          plot(x_apc_init,y_apc_init,'y--')
            hold(H,'off')
            p = get(gca,'position');
            h = axes('parent',gcf,'position',[p(1)+.45 p(2)+.6 0.3*p(3) 0.2*p(4)]);
            axis([0 60 0 12])
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
% %     scatter(1:n,ct(:,i))
% %     axis([0 60 0 0.5])
   F(idx_frame) = getframe(gcf);

% Record.pos(idx_frame,:,:) = pos;
% Record.Aforce(idx_frame,:) = alpha;
% Record.Bforce(idx_frame,:) = beta;
% Record.Lforce(idx_frame,:) = gamma;
    idx_frame = idx_frame + 1;
end

if mod(i,10000) == 0
    fprintf('progress: %d%% \n', i/N*100)
end


end
% sym = sqrt(centroid(1,:).^2 + centroid(2,:).^2);
% plot(1:N,sym)
% axis([1 N 0 0.1])

 toc
 %%
H = axes;
drawfig(H,pos);
hold(H,'on')
plot(x_vm,y_vm,'--')
hold(H,'off')
p = get(gca,'position');
h = axes('parent',gcf,'position',[p(1)+.45 p(2)+.6 0.3*p(3) 0.2*p(4)]);
axis([0 60 0 15])
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

%% 



%movie(gcf,F,1,10,[0 0 0 0])











