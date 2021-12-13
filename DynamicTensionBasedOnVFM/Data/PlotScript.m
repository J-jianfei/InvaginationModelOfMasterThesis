FileName = {'ApicalDynamicTension_ApicalConstrictionOnly','ApicalDynamicTension_LateralIncreaseNoBasalDepletion.mat','ApicalDynamicTension_BasalDepletionOnly.mat','ApicalDynamicTension_LateralIncreaseAndBasalDepletion.mat'};
i = 4; % 1: no lateral no basal; 2: lateral no basal; 3: no lateral but basal; 4:lateral and basal
load(FileName{i},'Record','x_vm','y_vm');
NumOfFrames = size(Record.pos,1);
NumOfCells = 0.5*size(Record.pos,3);

W = 0.5; % width of main axis
w = 0.1; % width of inset
% ParentAxis = axes('Position',[0.02 0.3 0.97 0.4]);
Col = 5;
for j = 1:Col+1
    if j == 1
        plotFrame = 1;
    else
        plotFrame = (j-1)/Col*NumOfFrames;
    end
    pos = squeeze(Record.pos(plotFrame,:,:));
    alpha = Record.Aforce(plotFrame,:);
    beta = Record.Bforce(plotFrame,:);
    gamma = Record.Lforce(plotFrame,:);
    process = plotFrame/NumOfFrames*100;
    
    posMainAxis = [-0.02+(j-1)*0.325*W-0.12 0.3 W W];
 %   posInset = [posMainAxis(1)+0.35*W, posMainAxis(2)+0.75*W, w, w];
   % h = axes('parent',gcf,'position',posMainAxis);
    h = axes('Position',posMainAxis);
    drawfig(h,pos);
%    text(-0.5+(j-1)*0.05*W ,0.0,"Process: "+num2str(process,'%.3f')+"%")
    hold(h,'on')
    plot(x_vm,y_vm,'--')
    hold(h,'off')
    axis("off")
%     Inset = axes('Position',posInset);
%     hold(Inset,"on")
%     plot(Inset,alpha,'r')
%     plot(Inset,beta,'b')
%     plot(Inset,gamma,'g')
%     axis([0 NumOfCells 0 15])
%     lgd = legend('apical','basal','lateral','location',[posInset(1)+.075 posInset(2)+.075 0.2*posInset(3) 0.2*posInset(4)]);
%     legend('boxoff')
%     title(lgd,'tensions')
%     xlabel('cell index')
%     hold(Inset,"off")
    
    
%     subplot("Position",positionInset)

  



end