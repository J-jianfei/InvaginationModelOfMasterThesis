clear 
clc
Filename = input("File name: ","s");
load(Filename)

numOfFig = 6;
frameInterval = size(Record.pos,1)/(numOfFig - 1);

W = 0.5; % width of main axis
w = 0.1; % width of inset

plotTension = 0;

for i = 1:numOfFig
    if i == 1
        frame = 1;
    else
        frame = (i-1)*frameInterval;
    end
    pos = squeeze(Record.pos(frame,:,:));
    alpha = squeeze(Record.Aforce(frame,:));
    beta = squeeze(Record.Bforce(frame,:));
    gamma = squeeze(Record.Lforce(frame,:));


    if plotTension == 0
        posMainAxis = [-0.02+(i-1)*0.325*W-0.12 0.3 W W];
        h = axes('Position',posMainAxis);
        drawfig(h,pos);
        hold(h,'on')
        plot(x_vm,y_vm,'--')
        hold(h,'off')
        axis("off")
    else
        W = 0.15;
        H = 0.5;
        posMainAxis = [0.02+(i-1)*1.1*W 0.3 W H];
        h = axes('Position',posMainAxis);
        plot(h,alpha,'r','LineWidth',2)
        hold(h,"on")
        plot(h,beta,'b','LineWidth',2)
        plot(h,gamma,'g','LineWidth',2)
        hold(h,"off")
        if i == 1
        legend("apical","basal","lateral")
        xlabel("cell index")
        ylabel("tension")
        else
            yticklabels("")
        end
        axis([0 80 0 12])
    end

end