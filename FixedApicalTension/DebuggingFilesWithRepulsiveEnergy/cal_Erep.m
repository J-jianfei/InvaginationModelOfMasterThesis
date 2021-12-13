function Erep = cal_Erep(pos,bas_l)
%calculate repulsive energy once basal vertices are inside another cell
% %% based on curvature

n = length(pos)/2;
E = zeros(1,n);
amp = 5;
L = 10^1;
% 
% posBas = pos(:,1:n);
% posPre = [posBas(:,end) posBas(:,1:end-1)];
% posNex = [posBas(:,2:n) posBas(:,1)];
% xDataCur = [posPre(1,:);posBas(1,:);posNex(1,:)];
% yDataCur = [posPre(2,:);posBas(2,:);posNex(2,:)];
% ta = sqrt((xDataCur(2,:)-xDataCur(1,:)).^2 + (yDataCur(2,:) - yDataCur(1,:)).^2);
% tb = sqrt((xDataCur(3,:)-xDataCur(2,:)).^2 + (yDataCur(3,:) - yDataCur(2,:)).^2);
% 
% M(1,1,:) = ones(1,n);
% M(1,2,:) = -ta;
% M(1,3,:) = ta.^2;
% M(2,1,:) = ones(1,n);
% M(2,2,:) = zeros(1,n);
% M(2,3,:) = zeros(1,n);
% M(3,1,:) = ones(1,n);
% M(3,2,:) = tb;
% M(3,3,:) = tb.^2;
% 
% a = zeros(3,n);
% b = zeros(3,n);
% for i = 1:n
%     a(:,i) = M(:,:,i)\xDataCur(:,i);
%     b(:,i) = M(:,:,i)\yDataCur(:,i);
% end
% 
% kappa = 2*(a(3,:).*b(2,:) - b(3,:).*a(2,:))./(a(2,:).^2 + b(2,:).^2).^(1.5);
% ProsVertex = find(abs(kappa)>=10);
%% based on basal length;
% 
if nargin == 2
  ProsVertex = find(bas_l<= 0.05);
if isempty(ProsVertex)
    Erep = 0;
else
    ProsVertex = [ProsVertex+1,ProsVertex-1];
    ProsVertex = unique(ProsVertex);
    
    for j = 1:length(ProsVertex) % find basal vertices inside the tissue
        % % % Serial code: very slow
        % %     for j = i-5:i+5
        % %         if j <= 0
        % %             J = j+n;
        % %         elseif j > n
        % %             J = j-n;
        % %         else
        % %             J = j;
        % %         end
        % %         k = J+1;
        % %         if k > n
        % %             k = k-n;
        % %         end
        % %         cell_vertices = pos(:,i);
        % %         A1 = polyarea([pos(1,i),pos(1,J),pos(1,k)],[pos(2,i) pos(2,J),pos(2,k)]);
        % %         A2 = polyarea([pos(1,i),pos(1,k),pos(1,k+n)],[pos(2,i) pos(2,k),pos(2,k+n)]);
        % %         A3 = polyarea([pos(1,i),pos(1,k+n),pos(1,J+n)],[pos(2,i) pos(2,k+n),pos(2,J+n)]);
        % %         A4 = polyarea([pos(1,i),pos(1,J+n),pos(1,J)],[pos(2,i) pos(2,J+n),pos(2,J)]);
        % %         Atot = A1+A2+A3+A4;
        % %         Acell = polyarea([pos(1,J),pos(1,k),pos(1,k+n),pos(1,J+n)],[pos(2,J),pos(2,k),pos(2,k+n),pos(2,J+n)]);
        % %         if abs(Atot - Acell)/Acell <= 10^-5
        % %             vec1 = [pos(:,i)-pos(:,J);0];
        % %             vec2 = [pos(:,k)-pos(:,J);0];
        % %             CrossProd = cross(vec1,vec2);
        % %             mag = norm(CrossProd);
        % %             d = mag/norm(vec2);
        % %             E(i) = E(i) + amp *(exp(L*d)-1);
        % %         end
        % %     end
        %% vectorized 
        %% calculate prospective vertices
        i = ProsVertex(j);
        CheckRange = i-5:i+5;
        CheckRange(CheckRange == i) =[];
        CheckRange(CheckRange == i-1) =[];
        CheckLength = length(CheckRange);
        CheckVertexBasal = [CheckRange;CheckRange+1];
        CheckVertexBasal(CheckVertexBasal>n) = CheckVertexBasal(CheckVertexBasal>n) - n;
        CheckVertexBasal(CheckVertexBasal<=0) = CheckVertexBasal(CheckVertexBasal<=0) + n;
        CheckVertexApical = [CheckVertexBasal(2,:);CheckVertexBasal(1,:)] + n;
        CheckVertex = [CheckVertexBasal;CheckVertexApical];
        xDataCheckPoint = pos(1,i)*ones(1,CheckLength);
        yDataCheckPoint = pos(2,i)*ones(1,CheckLength);
        xData1 = [xDataCheckPoint;pos(1,CheckVertex(1,:));pos(1,CheckVertex(2,:))];
        yData1 = [yDataCheckPoint;pos(2,CheckVertex(1,:));pos(2,CheckVertex(2,:))];
        A1 = polyarea(xData1,yData1);
        xData2 = [xDataCheckPoint;pos(1,CheckVertex(2,:));pos(1,CheckVertex(3,:))];
        yData2 = [yDataCheckPoint;pos(2,CheckVertex(2,:));pos(2,CheckVertex(3,:))];
        A2 = polyarea(xData2,yData2);
        xData3 = [xDataCheckPoint;pos(1,CheckVertex(3,:));pos(1,CheckVertex(4,:))];
        yData3 = [yDataCheckPoint;pos(2,CheckVertex(3,:));pos(2,CheckVertex(4,:))];
        A3 = polyarea(xData3,yData3);
        xData4 = [xDataCheckPoint;pos(1,CheckVertex(4,:));pos(1,CheckVertex(1,:))];
        yData4 = [yDataCheckPoint;pos(2,CheckVertex(4,:));pos(2,CheckVertex(1,:))];
        A4 = polyarea(xData4,yData4);
        Atot = A1 + A2 + A3 + A4;
        xDataCell = [pos(1,CheckVertex(1,:));pos(1,CheckVertex(2,:));pos(1,CheckVertex(3,:));pos(1,CheckVertex(4,:))];
        yDataCell = [pos(2,CheckVertex(1,:));pos(2,CheckVertex(2,:));pos(2,CheckVertex(3,:));pos(2,CheckVertex(4,:))];
        Acell = polyarea(xDataCell,yDataCell);
        Index = find(abs(Atot - Acell)./Acell <= 10^(-5));
        Len = length(Index);
        if Len ~= 0
        InsideWhichCell = CheckVertexBasal(1,Index);
        ApicalVertex = CheckVertexApical(2,Index);
        NextBasalVertex = CheckVertexBasal(2,Index);
        NextApicalVertex = CheckVertexApical(1,Index);
        %% calculate shortest deepness
        Vector1 = [pos(:,i).*ones(2,Len) - pos(:,InsideWhichCell);zeros(1,Len)];
        Vector2 = [pos(:,NextBasalVertex) - pos(:,InsideWhichCell);zeros(1,Len)];
        CrossProd1 = cross(Vector1,Vector2);
        mag = sqrt(CrossProd1(1,:).^2 + CrossProd1(2,:).^2 + CrossProd1(3,:).^2);
        Vec2Length = sqrt(Vector2(1,:).^2 + Vector2(2,:).^2 +Vector2(3,:).^2);
        d_fromBasal = mag./Vec2Length;
        % ********
        Vector3 = [pos(:,i).*ones(2,Len) - pos(:,InsideWhichCell);zeros(1,Len)];
        Vector4 = [pos(:,ApicalVertex) - pos(:,InsideWhichCell);zeros(1,Len)];
        CrossProd2 = cross(Vector3,Vector4);
        mag2 = sqrt(CrossProd2(1,:).^2 + CrossProd2(2,:).^2 + CrossProd2(3,:).^2);
        Vec4Length = sqrt(Vector4(1,:).^2 + Vector4(2,:).^2 +Vector4(3,:).^2);
        d_fromLateral_right = mag2./Vec4Length;
        % ********
        Vector5 = [pos(:,i).*ones(2,Len) - pos(:,NextBasalVertex);zeros(1,Len)];
        Vector6 = [pos(:,NextApicalVertex) - pos(:,NextBasalVertex);zeros(1,Len)];
        CrossProd3 = cross(Vector5,Vector6);
        mag3 = sqrt(CrossProd3(1,:).^2 + CrossProd3(2,:).^2 + CrossProd3(3,:).^2);
        Vec6Length = sqrt(Vector6(1,:).^2 + Vector6(2,:).^2 +Vector6(3,:).^2);
        d_fromLateral_left = mag3./Vec6Length;
        d = min([d_fromBasal;d_fromLateral_right;d_fromLateral_left],[],1);
        % energy
        E(i) = sum(amp*(exp(L*d)-1));
        else
            E(i) = 0;
        end
        
    end
    Erep = sum(E);

end

    
    
end

end

