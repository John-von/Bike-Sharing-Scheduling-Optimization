function [cost] = costFuction(curr_vc, a, b, s, L, dist, demands, cap, alpha, beta)
    [TD] = travel_distance(curr_vc, dist);
    [q] = violateLoad(curr_vc, demands, cap);
    [w] = violateTW(curr_vc, a, b, s, L, dist);
    cost = TD + alpha*q + beta*w;
end