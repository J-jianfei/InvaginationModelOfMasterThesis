clear all
clc
w = 0.1; % width of xyticks
W = 0.9; % width of parent axes;
Lines = 3;
Columns = 9;
figure
axes('position',[0.08 0.38 W 0.6*W])
xticks(0:w:0.9)
xticklabels({'','1.5','2.0','2.5', '3.0', '3.5','4.0','4.5','5.0','5.5'})
yticks(0:3*w:1)
yticklabels({'','0.5','1.0', '1.5', ''})
xlabel('Basal tension')
ylabel('Lateral tension')
title('Effect of Lateral and Basal Tension')
%axis('square')
p = get(gca,'position');

%% load data
namelist = dir('g*b*.mat');
l = length(namelist);
P = cell(1,l);

%% plot data
idx = 1;
for i = 1:Lines
    for j = 1:Columns
        h = axes('parent',gcf,'position',[p(1)+w*W*(j-1)+0.45*w p(2)+3*w*0.6*W*(i-1)+0.5*w w*W 3*w*0.6*W]);
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