clear 
clc
%% geometry-modified test version: no genetic regulation
tic

%% ***********  iteration parameters ************** %
idx_frame = 1;
num_frame = 400;
N =  10^5;
dt = 1.5*10^(-4);

%% default model coefficients (fixed)
n = 80; % number of cells for single layer embryo
R_y = 0.8; % raidus of yolk = radius of basal side of epithelium 
R_ap = 1.2; % radius of apical side of epithelium 
Ka = 10^4;   % cell area elasticity
Ky = 10^2; % yolk area elasticity
Epress = 10;
idx = 1:n;
x_bas = R_y*cos((idx - 1)*2*pi/n );
y_bas = R_y*sin((idx - 1)*2*pi/n);
x_apc = R_ap*cos((idx - 1)*2*pi/n);
y_apc = R_ap*sin((idx - 1)*2*pi/n);
pos = [x_bas x_apc;y_bas y_apc]; % initial position
yolk_area = polyarea(x_bas,y_bas);
Ay0 = yolk_area; % prefered yolk areads
cell_area = cal_cell_area(pos);
A0 = mean(cell_area); % prefered cell area
x_vm = [pos(1,n+1:end) pos(1,n+1)]; 
y_vm = [pos(2,n+1:end) pos(2,n+1)]; % get vitelline membrane position
r0 = mean(sqrt(x_vm(1:end-1).^2 + y_vm(1:end-1).^2)); % radius of membrane
Ad =10; % adhesion 
adsz = 15; % adhesion size
[ap_l,bas_l,lat_l]=cal_length(pos);
ap_l_0 =mean(ap_l);
d0 = 0.5*ap_l_0;
ita = 1; % drag factor
E = zeros(1,N); % Energy minimization records
% predefined regions of tissue
tissue = 1:n;
vm = 0.5*tissue(1) + 0.25*tissue(end); % predefined ventral midline
center = [floor(vm) ceil(vm)]; % cell index at the midline
distance = tissue - vm; 
distance(distance>0.5*n) = distance(distance>0.5*n) - n;




%%  fixed Tension
% mu = 0.02;
% alpha0 = 9;
% beta0 = 5.0;
% gamma0 = 1.0;
% baseline = 1;
% 
% alpha = alpha0*exp(-mu*distance.^2)+baseline;
% beta = beta0 * ones(1,n);
% gamma = gamma0 * ones(1,n);

%% genetic regulation parameters

r1 = 0.75; % relocalization coeff from basal to apical
r2 = 0.15;
ut = 0.8; % twist breakdown rate;
um = 1; % myosin breakdown rate;
pt = 1; % twist production rate;
pm = 0.0; % constant myosin production rate;
ka = 1; % saturation of twi
lm = 0.03; % basal myosin loss rate on mesoderm
pm_e = 0.25; % basal myosin production rate
um_e = 0.1; % basal myosin breakdown rate 
mu = 0.02; % std of Gaussian term
kh =1;
pm_l = 0.75    ; % lateral constant myosin production
Th = 0.3; % the value when half-threshold of overall production term is reached 
Qm = 0.10; % the value when half-threshold of mechanical term is reached
q = 3; % power in hill type mechanical response
l = 8; % hill term of mechanical response
pth = 1.5;
a = 0.5; % contribution from genetic term
b = 0.5; % contribution from mechnical term
fb = 1; % basal myosin concentration to tension factor
fa = 1; % apical myosin concentration to tension factor
fg = 1; % lateral myosin concentration to tension factor
%*********** initial conditons and constraints *******%
cm_apc = zeros(n,N);
cm_apc(:,1) = 0.0;
ct = zeros(n,N);
cm_bas = zeros(n,N);
cm_bas(:,1) = 4;
cm_lat =zeros(n,N);
alpha_min = 0.5;
beta_min = 0.1;
gamma_min = 0.1;

%******** mutant experiment **********%
mutant = 0; % 0:normal; 1:global elliptic compression; 2:ectopic compression
if mutant == 1
    load('pos_ellipse20.mat','pos');
    x_vm = pos(1,n+1:end);
    y_vm = pos(2,n+1:end);
    minor_axis = 0.8*r0;
    major_axis = 1.2*r0;
    [ap_l,bas_l,lat_l] = cal_length(pos);
elseif mutant == 2
    load('pos_ecto.mat','pos');
    [ap_l,bas_l,lat_l] = cal_length(pos);
end



Delay = 0.0*N; % consider delayed lateral increase
for i = 1:N
    %% concentration-to-tension

    alpha = fa*cm_apc(:,i)';
    beta = fb*cm_bas(:,i)';
    gamma = fg*cm_lat(:,i)';
    alpha(alpha<alpha_min) = alpha_min;
    beta(beta<beta_min) = beta_min;
    gamma(gamma<gamma_min) = gamma_min;

    %% energy minimization
    if mutant == 1
        [Egrad,Etot] = cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,Ay0,A0,center,d0,Ad,adsz,Epress,minor_axis,major_axis);
    else
        [Egrad,Etot] = cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,Ay0,A0,center,d0,Ad,adsz,Epress,r0);
    end
    net_force = - Egrad;
    E(i) = Etot;
    disp_f = dt/ita*net_force;
    pos = pos + disp_f;
    [ap_l,bas_l,lat_l]= cal_length(pos);
    %% genetic regulation equations
    mech_term = 1./((Qm./((ap_l' - ap_l_0)/ap_l_0)).^l + 1);
    gen_term = exp(-mu*distance'.^2);
    prod = b*mech_term + a*gen_term;
    HillTerm = pth./(((Th./prod).^q)+1);
    dct = dt*(pt*ct(:,i).^2./(1+ka*ct(:,i).^2) - ut*ct(:,i) + HillTerm);
    dcm_apc = dt*(r1*ct(:,i).^2- um*cm_apc(:,i)+pm);
    dcm_bas = dt * ( - lm*ct(:,i).^2 + pm_e*cm_bas(:,i) - um_e*cm_bas(:,i).^2);
    %  dcm_lat = dt * ( r2*ct(:,i).^2 - um*cm_lat(:,i) + pm_l);
    % dcm_lat = dt * ( r2*ct(:,i).^2 - um*cm_lat(:,i).^2 + pm_l_log*cm_lat(:,i));
    dcm_lat = dt * ( r2*ct(:,i).^2 - um*cm_lat(:,i) + pm_l);

    ct(:,i+1) = ct(:,i) + dct;
    ct(ct(:,i+1)<=0,i+1) = 0 ;
    cm_apc(:,i+1) = cm_apc(:,i)+dcm_apc;
    cm_bas(:,i+1) = cm_bas(:,i) + dcm_bas;
    cm_lat(:,i+1) = cm_lat(:,i) + dcm_lat;




%% movie 
if mod(i,num_frame)==1    
    Record.pos(idx_frame,:,:) = pos;
    Record.Aforce(idx_frame,:) = alpha;
    Record.Bforce(idx_frame,:) = beta;
    Record.Lforce(idx_frame,:) = gamma;
    idx_frame = idx_frame + 1;
end


    

 if mod(i,10000) == 0
     clc
     fprintf('progress: %d%% \n', i/N*100)
     toc
 end

end
%%
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



%movie(gcf,F,1,10,[0 0 0 0])




