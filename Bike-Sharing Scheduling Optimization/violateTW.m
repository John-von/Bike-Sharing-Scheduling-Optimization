function [w] = violateTW(curr_vc, a, b, s, L, dist)
    NV = size(curr_vc, 1);
    w = 0;
    bsv = begin_s_v(curr_vc, a, s, dist);
    for i = 1:NV
        route = curr_vc{i};
        bs = bsv{i};
        l_bs = length(bsv{i});
        for j = 1:l_bs-1
            if bs(j) > b(route(j))
                w = w + bs(j) - b(route(j));
            end
        end
        if bs(end) > L
            w = w + bs(end) - L;
        end
    end
end