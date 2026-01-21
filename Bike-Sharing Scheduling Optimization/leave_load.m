function [Ld, initial_load] = leave_load(route, demands)
    % Ld: 路线上需要的最大装载量
    % initial_load: 从调度中心出发时的初始装载
    
    n = length(route);
    if n == 0
        Ld = 0;
        initial_load = 0;
        return;
    end
    
    current_load = 0;  % 当前车上的自行车数
    max_load = 0;      % 沿途最大装载量
    
    for i = 1:n
        station = route(i);
        if demands(station) > 0
            % 站点需要车（卸货）
            current_load = current_load - demands(station);
        end
        max_load = max(max_load, current_load);
    end
    
    Ld = max_load;
    initial_load = max_load;  % 初始需要装载的量
end