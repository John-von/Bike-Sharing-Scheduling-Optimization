function [q] = violateLoad(curr_vc, demands, cap)
    NV = size(curr_vc, 1);
    q = 0;
    for i = 1:NV
        route = curr_vc{i};
        Ld = leave_load(route, demands);
        if Ld > cap
            q = q + Ld - cap;
        end
    end
end