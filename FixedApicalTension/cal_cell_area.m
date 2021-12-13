function cell_area = cal_cell_area(pos)
% this function aims to calculate cell area at each iteration process.
% input pos is the position of vertices, indexed in the same convention as
% pervious
num_v = length(pos); % number of vertices
num_c = 0.5*num_v;  % number of cells
ap_v = pos(:,num_c+1:end);
bas_v = pos(:,1:num_c);
x = zeros(4,num_c);
y = zeros(4,num_c);

x(:,1:end-1) = [ap_v(1,1:end-1);ap_v(1,2:end);bas_v(1,2:end);bas_v(1,1:end-1)];
x(:,end) = [ap_v(1,end);ap_v(1,1);bas_v(1,1);bas_v(1,end)];

y(:,1:end-1) = [ap_v(2,1:end-1);ap_v(2,2:end);bas_v(2,2:end);bas_v(2,1:end-1)];
y(:,end) = [ap_v(2,end);ap_v(2,1);bas_v(2,1);bas_v(2,end)];

cell_area = polyarea(x,y);

% for i = 1:num_c
%     if i ~= num_c
%     x = [pos(1,i) pos(1,i+num_c) pos(1,i+1+num_c) pos(1,i+1)];
%     y = [pos(2,i) pos(2,i+num_c) pos(2,i+1+num_c) pos(2,i+1)];
%     cell_area(i) = polyarea(x,y);
%     else
%     x = [pos(1,i) pos(1,i+num_c) pos(1,i+1) pos(1,1)];
%     y = [pos(2,i) pos(2,i+num_c) pos(2,i+1) pos(2,1)];
%     cell_area(i) = polyarea(x,y);
%     end
% 
% end
end

