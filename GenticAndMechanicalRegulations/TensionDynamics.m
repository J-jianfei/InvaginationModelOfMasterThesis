n = 80;
% iteration variables
N = 10^5;
dt = 1.5*10^-4;
% predefined regions of tissue
tissue = 1:n;
vm = 0.5*tissue(1) + 0.25*tissue(end); % predefined ventral midline
distance = tissue - vm; 
distance(distance>0.5*n) = distance(distance>0.5*n) - n;

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
pm_l = 1    ; % lateral constant myosin production
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

Delay = 0;
for i = 1:N


    gen_term = exp(-mu*distance'.^2);
    dct = dt*(pt*ct(:,i).^2./(1+ka*ct(:,i).^2) - ut*ct(:,i) + pth*gen_term);
    dcm_apc = dt*(r1*ct(:,i).^2- um*cm_apc(:,i)+pm);
    dcm_bas = dt * ( - lm*ct(:,i).^2 + pm_e*cm_bas(:,i) - um_e*cm_bas(:,i).^2);
    iDelay = i - Delay;
    if iDelay >= 1
    dcm_lat = dt * ( r2*ct(:,i).^2 - um*cm_lat(:,i) + pm_l);
    else
        dcm_lat = 0;
    end

    ct(:,i+1) = ct(:,i) + dct;
    ct(ct(:,i+1)<=0,i+1) = 0 ;
    cm_apc(:,i+1) = cm_apc(:,i)+dcm_apc;
    cm_bas(:,i+1) = cm_bas(:,i) + dcm_bas;
    cm_lat(:,i+1) = cm_lat(:,i) + dcm_lat;
end

%% apical tension dynamics: spatial resolved
subplot(121)
plot(cm_apc(:,1))
hold on
plot(cm_apc(:,N/10))
plot(cm_apc(:,N/4))
plot(cm_apc(:,N/2))
plot(cm_apc(:,end))
hold off
lgd = legend("t = 0","t = 0.1N","t = 0.25N","t = 0.5N","t = N");
xlabel("cell index")
ylabel("tension")

%% tension dynamics at ventral midline
subplot(122)
plot(cm_apc(20,:),'r')
hold on
plot(cm_bas(20,:),'b')
plot(cm_lat(20,:),'g')
hold off
axis([0 10^5 0 10])
legend("apical force","basal force", "lateral force");
xlabel("iteration process")
ylabel("tension")
