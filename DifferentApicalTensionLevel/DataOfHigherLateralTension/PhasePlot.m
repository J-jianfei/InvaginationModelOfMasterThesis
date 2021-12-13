clear all
clc
w = 0.2; % width of xyticks
W = 0.8; % width of parent axes;
Lines = 4;
Columns = 4;
figure
axes('position',[0.08 0.08 W W])
xticks(0:w:1)
xticklabels({'','1.5','2.5','4.5', '6.5', ''})
yticks(0:w:1)
yticklabels({'','1.5','2.5','4.5', '6.5', ''})
xlabel('Apical tension')
ylabel('Basal tension')
title('Effect of Increased Apical/Basal Tension')
%axis('square')
p = get(gca,'position');

%% load data
namelist = dir('a*b*g2.5.mat');
l = length(namelist);
P = cell(1,l);

%% plot data
idx = 1;
for i = 1:Lines
    for j = 1:Columns
        h = axes('parent',gcf,'position',[p(1)+w*W*(j-1/2)-0.04 p(2)+w*W*(i-0.5)-0.05 1.5*w*W 1.5*w*W]);
        p1 = get(gca,'position');
        filename = namelist(idx).name;
        P{idx} = load(filename,'pos');
        drawfig(gca,P{idx}.pos);
        % box off
        axis off
        xticks('');
        yticks('');
        idx = idx + 1;
    end
end