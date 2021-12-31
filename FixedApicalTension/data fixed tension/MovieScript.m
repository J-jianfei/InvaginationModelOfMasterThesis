clear all
clc
% load('Record_3.8.mat')
% load('Record_3.75.mat')
% load('Record_normal_3.75.mat')
% load('Record_normal_3.8.mat')
% load('Record_local_compress_3.8.mat')
% load('Record_local_compress_3.75_N.mat')
% load('b0.5g1.0rep.mat')
% load('Demo.mat')
% load('bug2.mat');
load('g1.0b2.5.mat')
NumofFrames = size(Record.pos,1);
NumofVertices = size(Record.pos,3);

% NumofFrames = 100;
for i = 1:NumofFrames
    pos = squeeze(Record.pos(i,:,:));
    if i == 1
        x_vm = pos(1,NumofVertices/2+1:end);
        y_vm = pos(2,NumofVertices/2+1:end);
    end
    alpha = squeeze(Record.Aforce(i,:));
    beta = squeeze(Record.Bforce(i,:));
    gamma = squeeze(Record.Lforce(i,:));
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
    F(i) =getframe(gcf);
end
%%
v = VideoWriter('GeoCorrectedDemo');
v.FrameRate = 50;
open(v)
writeVideo(v,F)
close(v)