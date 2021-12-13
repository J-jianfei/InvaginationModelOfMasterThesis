function [E_grad,Etot] = cal_E_grad(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0,r1) 
%calculate energy gradient
%


    n = length(pos);
    step_length = 10^(-9); % delta
%% nested loop
%    num_v = length(pos); % number of vertices
%     E_grad = zeros(2,num_v);
%     
%     if nargin == 17
%         E0 = cal_Etot(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0);
%         Etot = E0;
%         for j = 1:2
%             parfor i = 1:num_v
%                 pos_temp = pos;
%                 pos_temp(j,i) = pos(j,i)+ step_length; % x'= x + dx;
%                 [ap_l_temp,bas_l_temp,lat_l_temp] = cal_length(pos_temp);
%                 E = cal_Etot(ap_l_temp,bas_l_temp,lat_l_temp,pos_temp,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0);
%                 E_grad(j,i) = (E - E0)/step_length;
%             end
%         end
%     elseif nargin == 18
%         E0 = cal_Etot(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0,r1);
%         for j = 1:2
%             parfor i = 1:num_v
%                 pos_temp = pos;
%                 pos_temp(j,i) = pos(j,i)+ step_length; % x'= x + dx;
%                 [ap_l_temp,bas_l_temp,lat_l_temp] = cal_length(pos_temp);
%                 E = cal_Etot(ap_l_temp,bas_l_temp,lat_l_temp,pos_temp,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0,r1);
%                 E_grad(j,i) = (E - E0)/step_length;
%             end
%         end
%     end


%% vectorized

NumOfInput =nargin;

xData = pos(1,:);
yData = pos(2,:);
xData = repmat(xData,[n,1]);
yData = repmat(yData,[n,1]);
dl = step_length*eye(n);
xData_plus_dx = xData + dl ;
yData_plus_dy = yData + dl;

if NumOfInput == 18
    Etot = cal_Etot(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0,r1);
    for i = 1:n
        pos_temp_x = [xData_plus_dx(i,:);yData(i,:)];
        pos_temp_y = [xData(i,:);yData_plus_dy(i,:)];
        [ap_l_temp_x,bas_l_temp_x,lat_l_temp_x] = cal_length(pos_temp_x);
        [ap_l_temp_y,bas_l_temp_y,lat_l_temp_y] = cal_length(pos_temp_y);
        E_temp_x = cal_Etot(ap_l_temp_x,bas_l_temp_x,lat_l_temp_x,pos_temp_x,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0,r1);
        E_temp_y = cal_Etot(ap_l_temp_y,bas_l_temp_y,lat_l_temp_y,pos_temp_y,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0,r1);
        E_grad(:,i) = [(E_temp_x - Etot)/step_length;(E_temp_y - Etot)/step_length];
    end  
elseif NumOfInput == 17
    Etot = cal_Etot(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0);
    for i = 1:n
        pos_temp_x = [xData_plus_dx(i,:);yData(i,:)];
        pos_temp_y = [xData(i,:);yData_plus_dy(i,:)];
        [ap_l_temp_x,bas_l_temp_x,lat_l_temp_x] = cal_length(pos_temp_x);
        [ap_l_temp_y,bas_l_temp_y,lat_l_temp_y] = cal_length(pos_temp_y);
        E_temp_x = cal_Etot(ap_l_temp_x,bas_l_temp_x,lat_l_temp_x,pos_temp_x,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0);
        E_temp_y = cal_Etot(ap_l_temp_y,bas_l_temp_y,lat_l_temp_y,pos_temp_y,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,adsz,press,r0);
        E_grad(:,i) = [(E_temp_x - Etot)/step_length;(E_temp_y - Etot)/step_length];
    end    
end
%     

    
end

