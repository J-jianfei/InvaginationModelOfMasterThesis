clear all
clc
load('demo.mat')
[Egrad,Etot] =  cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,0.5*ap_l_0,adhesion,adsz,press,r0);