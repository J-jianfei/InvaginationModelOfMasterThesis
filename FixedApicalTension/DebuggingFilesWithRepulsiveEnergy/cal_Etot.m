function Etot = cal_Etot(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0,center,ref_l,adh,ad_sz,press,r0,r1)
% calculate total energy
 NumOfInput =nargin;
    if NumOfInput == 17
        Eext = cal_Eext(press,pos,r0);
    elseif NumOfInput == 18
        Eext = cal_Eext(press,pos,r0,r1);
    end
    E = cal_E(ap_l,bas_l,lat_l,pos,alpha,beta,gamma,Ka,Ky,ya0,ca0);
    Ead = cal_Ead(pos,center,ref_l,adh,ad_sz);
    Erep = cal_Erep(pos,bas_l);
    Etot = Ead + E + Eext + Erep;

% if isnan(Etot)
%     if isnan(E)
%         fprintf('cell and yolk energy is NaN, check parameters\n')
%     end
%     if isnan(Eext)
%         fprintf('external vertices energy is NaN, check parameters\n')
%     end
%     if isnan(Ead)
%         fprintf('adhesion energy is NaN, check parameters\n')
%     end
%     pause
% end
end

