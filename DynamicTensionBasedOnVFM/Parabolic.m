function Tens = Parabolic(sz,val_ref,type)
% returns the parabolic tension profile 
% input: tissue size and maximal value of the tension,type: concave or
% convex parabolic relation
    tissue = -0.5*(sz-1):0.5*(sz-1);
if type == 1
    Tens = - val_ref/(0.5*(sz-1))^2 * tissue.^2 + val_ref + 1;
elseif type == 2
    Tens = (val_ref - 0.5)/(0.5*(sz-1))^2 * tissue.^2 + 0.5; 
end
end

