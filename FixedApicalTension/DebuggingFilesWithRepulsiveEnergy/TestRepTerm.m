clear all
% clc
load('b0.5g1.0.mat')

i = 30;
pos = squeeze(Record.pos(i,:,:));
alpha = squeeze(Record.Aforce(i,:));
beta = squeeze(Record.Bforce(i,:));
gamma = squeeze(Record.Lforce(i,:));
tic
for ii = 1:100
[ap_l,bas_l,lat_l]=cal_length(pos);
net_force = - cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,r0);
disp_f = dt/ita*net_force;
pos = pos + disp_f;

end
toc