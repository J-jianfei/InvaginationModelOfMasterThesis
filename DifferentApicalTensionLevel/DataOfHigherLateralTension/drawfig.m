function drawfig(ax,vertex_position)
% this function aims to draw the vertex configuration with input:
% vertex_position, which is a 2*2n matrix where n is the number of cells
% indices from 1 to n are basal vertices; indices from n+1 to 2n are apical vertices
% num_v = length(vertex_position); % number of vertices
% num_c = num_v/2; % number of cells
% xycoord = [vertex_position(1,:);vertex_position(2,:)]';
% A = zeros(num_v,num_v); % difine adjecent martix, connected vertices(i,j) with A(i,j)=A(j,i)=1; else is zero
% 
% for i = 1:num_v
%     if i < num_c
%         A(i,i+1) = 1;
%         A(i+1,i) = 1;
%         A(i,i+num_c) = 1;
%         A(num_c+i,i) = 1;
%     elseif i == num_c
%         A(1,num_c) = 1;
%         A(num_c,1) = 1;
%         A(i,i+num_c) = 1;
%         A(num_c+i,i) = 1;
%     elseif i > num_c && i < num_v
%         A(i,i+1) = 1;
%         A(i+1,i) = 1;
%     elseif i == num_v
%         A(num_c+1,num_v) = 1;
%         A(num_v,num_c+1) = 1;         
%     end 
% end
% 
% gplot(A,xycoord,'-')

n = length(vertex_position); % number of vertices
N = n/2; % number of cells

xdata = vertex_position(1,:);
ydata = vertex_position(2,:);

s = [1:N N+1:n 1:N];
t = [2:N 1 N+2:n N+1 N+1:n];

G = graph(s,t);

P = plot(ax,G,'-','XData',xdata,'YData',ydata,'Marker','none','LineWidth',2);
labelnode(P,1:n,'')

% xlim([-1.5 2.5])
% ylim([-1.5 2.5])
xlim([-2 2])
ylim([-2 2])
axis square

end

