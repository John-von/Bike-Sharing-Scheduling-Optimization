function route_table = generate_route_table(VC, demands, dist, cap, a, b, L, s, ...
                                              fixed_cost, var_cost, service_cost, ...
                                              hourly_cost, speed, actual_idx)
    if nargin < 14
        actual_idx = [];
    end
    
    NV = length(VC);
    RouteID = (1:NV)';
    Stations = cell(NV, 1);
    NumStations = zeros(NV, 1);
    TotalDemand = zeros(NV, 1);
    RouteDistance = zeros(NV, 1);
    RouteCost = zeros(NV, 1);
    LoadRate = zeros(NV, 1);
    ViolatesTW = zeros(NV, 1);
    
    for i = 1:NV
        route = VC{i};
        if isempty(actual_idx)
            actual_route = route;
        else
            actual_route = actual_idx(route);
        end
        
        % ========== 修改部分：生成带箭头的站点序列 ==========
        if isempty(actual_route)
            Stations{i} = 'Depot → Depot';
        else
            % 构建站点序列字符串：Depot → 1 → 5 → 8 → 12 → Depot
            station_str = 'Depot';
            for j = 1:length(actual_route)
                station_str = sprintf('%s → %d', station_str, actual_route(j));
            end
            station_str = sprintf('%s → Depot', station_str);
            Stations{i} = station_str;
        end
        % ====================================================
        
        NumStations(i) = length(route);
        TotalDemand(i) = sum(demands(route));
        LoadRate(i) = TotalDemand(i) / cap;
        
        % 计算路线距离
        route_dist = 0;
        if ~isempty(route)
            route_dist = dist(1, route(1)+1);
            for j = 2:length(route)
                route_dist = route_dist + dist(route(j-1)+1, route(j)+1);
            end
            route_dist = route_dist + dist(route(end)+1, 1);
        end
        RouteDistance(i) = route_dist;
        
        % 计算路线成本
        travel_time = route_dist / speed;
        RouteCost(i) = fixed_cost + var_cost * route_dist + ...
                       hourly_cost * travel_time + service_cost * NumStations(i);
        
        % 检查时间窗违反
        if ~isempty(route)
            [bs, back] = begin_s(route, a, s, dist);
            b_route = b(route);
            violate_tw = any(bs(:) > b_route(:)) || (back > L);
            ViolatesTW(i) = double(violate_tw);
        end
    end
    
    route_table = table(RouteID, Stations, NumStations, TotalDemand, ...
                        RouteDistance, RouteCost, LoadRate * 100, ViolatesTW, ...
                        'VariableNames', {'RouteID', 'Stations', 'NumStations', ...
                        'TotalDemand', 'RouteDistance_km', 'RouteCost_CNY', ...
                        'LoadRate_percent', 'ViolatesTW'});
end