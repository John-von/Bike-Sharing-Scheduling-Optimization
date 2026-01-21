function [bs, back] = begin_s(route, a, s, dist)
    n = length(route);
    bs = zeros(1, n);
    if n > 0
        bs(1) = max(a(route(1)), dist(1, route(1)+1));
        for i = 2:n
            bs(i) = max(a(route(i)), bs(i-1) + s(route(i-1)) + dist(route(i-1)+1, route(i)+1));
        end
        back = bs(end) + s(route(end)) + dist(route(end)+1, 1);
    else
        back = 0;
    end
end