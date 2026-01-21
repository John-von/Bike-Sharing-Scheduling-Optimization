function load_rates = calculate_load_rates(VC, demands, cap)
    NV = length(VC);
    load_rates = zeros(NV, 1);
    for i = 1:NV
        route = VC{i};
        if ~isempty(route)
            total_demand = sum(demands(route));
            load_rates(i) = total_demand / cap;
        end
    end
end