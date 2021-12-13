function T = ParaTime(time,t0,Tmax)
% parabolic time dependence of tensions
% input: time point it reaches its maximum t0, and maximum value Tmax,and
% current time (could be a vector)
% Note: Tmax here is the maximum in time, not in space.
T = -Tmax/t0^2*(time - t0).^2 + Tmax ;
end

