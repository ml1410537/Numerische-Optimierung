function out = solver_shiftcontrol(tvec, uvec, dt)
%SOLVER_SHIFTCONTROL Shift control trajectory by sampling time dt
    out = interp1([tvec tvec(end)+dt], [uvec uvec(:,end)]', tvec+dt);
    if size(uvec, 1) > 1
        out = out';
    end
end