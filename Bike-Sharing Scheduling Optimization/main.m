% =========================================================================
% Bike-Sharing Scheduling Optimization - Academic Figure System
% =========================================================================
clear; clc; close all; tic

fprintf('=================================================================\n');
fprintf('Bike-Sharing Hybrid Scheduling - Complete Academic System\n');
fprintf('=================================================================\n\n');

% =========================================================================
% Visual Configuration
% =========================================================================
FONT_NAME = 'Times New Roman';
FONT_SIZE_TITLE = 11;
FONT_SIZE_SUBPLOT = 10;
FONT_SIZE_LABEL = 9;
FONT_SIZE_LEGEND = 8;
FONT_SIZE_TICK = 8;

LINE_WIDTH_MAIN = 1.5;
LINE_WIDTH_GRID = 0.5;
LINE_WIDTH_AXIS = 1.0;
MARKER_SIZE = 7;

COLOR_SCHEME = struct();
COLOR_SCHEME.blue = [0, 0.447, 0.741];
COLOR_SCHEME.red = [0.850, 0.325, 0.098];
COLOR_SCHEME.green = [0.466, 0.674, 0.188];
COLOR_SCHEME.orange = [0.929, 0.694, 0.125];
COLOR_SCHEME.purple = [0.494, 0.184, 0.556];
COLOR_SCHEME.depot = [0.2, 0.2, 0.2];
COLOR_SCHEME.grid = [0.85, 0.85, 0.85];

set(0, 'DefaultAxesFontName', FONT_NAME);
set(0, 'DefaultAxesFontSize', FONT_SIZE_TICK);
set(0, 'DefaultTextFontName', FONT_NAME);
set(0, 'DefaultAxesLineWidth', LINE_WIDTH_AXIS);
set(0, 'DefaultAxesGridAlpha', 0.15);
set(0, 'DefaultAxesBox', 'on');

% =========================================================================
% Create Output Directories
% =========================================================================
output_dir = 'academic_outputs';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n\n', output_dir);
end

figures_dir = fullfile(output_dir, 'figures');
tables_dir = fullfile(output_dir, 'tables');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end

% =========================================================================
% Data Loading
% =========================================================================
fprintf('Loading data...\n');
c101 = importdata('c101.txt');
cap = 50;  % Vehicle capacity
E = c101(1,5);
L = c101(1,6);
vertexs = c101(:,2:3);
customer = vertexs(2:end,:);
cusnum = size(customer,1);
v_num = 10;  % Number of vehicles

demands = c101(2:end,4);
a = c101(2:end,5);  % Time window start
b = c101(2:end,6);  % Time window end
s = c101(2:end,7);  % Service time

h = pdist(vertexs);
dist = squareform(h);

% Cost parameters
fixed_cost_per_vehicle = 200;
variable_cost_per_km = 2.5;
service_cost_per_station = 8;
driver_hourly_cost = 50;
avg_speed = 30;  % km/h

fprintf('Data loaded: %d stations\n\n', cusnum);

%% =========================================================================
% Figure 1: Station Distribution Map
% =========================================================================
fprintf('Generating Figure 1: Station Distribution...\n');
figure('Position', [50, 50, 1000, 750], 'Color', 'w', 'PaperPositionMode', 'auto');
hold on; box on;

% Voronoi diagram background
[vx, vy] = voronoi(customer(:,1), customer(:,2));
plot(vx, vy, '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);

% Connecting lines to depot (sparse to avoid clutter)
for i = 1:5:cusnum
    plot([vertexs(1,1), customer(i,1)], [vertexs(1,2), customer(i,2)], ...
         ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
end

% Stations colored by distance from depot
distances = dist(1, 2:end)';
scatter(customer(:,1), customer(:,2), 250, distances, 'o', 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.8);
colormap(jet);
cb = colorbar;
ylabel(cb, 'Distance from depot (km)', 'FontSize', FONT_SIZE_LEGEND);

% Depot marker
plot(vertexs(1,1), vertexs(1,2), 'p', 'MarkerSize', 25, ...
     'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 2.5);

% Station ID labels
for i = 1:cusnum
    text(customer(i,1), customer(i,2), sprintf('%d', i), ...
         'FontSize', 8, 'Color', 'w', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'EdgeColor', 'none');
end

% Distance circles
theta_circle = linspace(0, 2*pi, 100);
distances_circles = [20, 40, 60];
circle_colors = [0.7 0.7 0.7; 0.6 0.6 0.6; 0.5 0.5 0.5];
for i = 1:length(distances_circles)
    d = distances_circles(i);
    x_circle = vertexs(1,1) + d * cos(theta_circle);
    y_circle = vertexs(1,2) + d * sin(theta_circle);
    plot(x_circle, y_circle, '--', 'Color', circle_colors(i,:), 'LineWidth', 1.2);
    text(vertexs(1,1) + d*0.707, vertexs(1,2) + d*0.707, sprintf('%d km', d), ...
         'FontSize', FONT_SIZE_LEGEND, 'Color', circle_colors(i,:), ...
         'FontWeight', 'bold', 'BackgroundColor', 'w');
end

xlabel('X coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Y coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridLineStyle', '-', 'GridAlpha', 0.3);
axis equal tight;

% Information box
area_size = (max(customer(:,1))-min(customer(:,1))) * (max(customer(:,2))-min(customer(:,2)));
avg_distance = mean(distances);
max_distance = max(distances);
info_text = sprintf(['Total Stations: %d\n' ...
                     'Depot: (%.1f, %.1f)\n' ...
                     'Service Area: %.0f km²\n' ...
                     'Avg Distance: %.1f km\n' ...
                     'Max Distance: %.1f km'], ...
                    cusnum, vertexs(1,1), vertexs(1,2), area_size, ...
                    avg_distance, max_distance);
text(0.02, 0.98, info_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', FONT_SIZE_LEGEND, ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, ...
     'Margin', 5);
xlim([115,122]);ylim([21,27]);
fig1_path = fullfile(figures_dir, 'Fig1_Station_Distribution');
savefig(gcf, [fig1_path '.fig']);
fprintf('  Saved: Fig1_Station_Distribution.fig\n');

%% =========================================================================
% Step 1: Baseline Algorithm Optimization (Operator-only dispatch)
% =========================================================================
fprintf('\nStep 1: Baseline Algorithm Optimization (Operator-only)\n');
fprintf('-----------------------------------------------------------------\n');

alpha = 50;
beta = 100;
NIND = 100;  % Population size
MAXGEN = 50;  % Max generations
Pc = 0.9;  % Crossover probability
Pm = 0.05;  % Mutation probability
GGAP = 0.9;  % Generation gap

N = cusnum + v_num - 1;
init_vc = init(cusnum, a, demands, cap);
Chrom = InitPopCW(NIND, N, cusnum, init_vc);

% Record iteration history
baseline_history = struct();
baseline_history.generation = zeros(MAXGEN, 1);
baseline_history.best_obj = zeros(MAXGEN, 1);
baseline_history.avg_obj = zeros(MAXGEN, 1);
baseline_history.worst_obj = zeros(MAXGEN, 1);

fprintf('Baseline optimization in progress...');
gen = 1;
while gen <= MAXGEN
    ObjV = calObj_baseline(Chrom, cusnum, cap, demands, a, b, L, s, dist, alpha, beta);
    
    baseline_history.generation(gen) = gen;
    baseline_history.best_obj(gen) = min(ObjV);
    baseline_history.avg_obj(gen) = mean(ObjV);
    baseline_history.worst_obj(gen) = max(ObjV);
    
    FitnV = Fitness(ObjV);
    SelCh = Select(Chrom, FitnV, GGAP);
    SelCh = Recombin(SelCh, Pc);
    SelCh = Mutate(SelCh, Pm);
    SelCh = LocalSearch_baseline(SelCh, cusnum, cap, demands, a, b, L, s, dist, alpha, beta);
    Chrom = Reins(Chrom, SelCh, ObjV);
    Chrom = deal_Repeat(Chrom);
    
    if mod(gen, 10) == 0
        fprintf('.');
    end
    gen = gen + 1;
end
fprintf('Complete\n');

ObjV = calObj_baseline(Chrom, cusnum, cap, demands, a, b, L, s, dist, alpha, beta);
[~, minInd] = min(ObjV);
[baseline_VC, baseline_NV, baseline_TD, baseline_violate_num, baseline_violate_cus] = ...
    decode_baseline(Chrom(minInd(1),:), cusnum, cap, demands, a, b, L, s, dist);

baseline_time_hours = baseline_TD / avg_speed;
baseline_vehicle_cost = baseline_NV * fixed_cost_per_vehicle;
baseline_distance_cost = baseline_TD * variable_cost_per_km;
baseline_time_cost = baseline_time_hours * driver_hourly_cost;
baseline_service_cost = cusnum * service_cost_per_station;
baseline_total_cost = baseline_vehicle_cost + baseline_distance_cost + ...
                      baseline_time_cost + baseline_service_cost;

fprintf('\nBaseline Results:\n');
fprintf('  Vehicles: %d\n', baseline_NV);
fprintf('  Total Distance: %.2f km\n', baseline_TD);
fprintf('  Violated Routes: %d\n', baseline_violate_num);
fprintf('  Violated Stations: %d\n', baseline_violate_cus);
fprintf('  Total Cost: ¥%.2f\n', baseline_total_cost);

%% =========================================================================
% Figure 2: Baseline Algorithm Convergence Curve
% =========================================================================
fprintf('\nGenerating Figure 2: Baseline Convergence...\n');
figure('Position', [100, 50, 1000, 750], 'Color', 'w', 'PaperPositionMode', 'auto');
hold on; box on;

% Fill area between best and worst
fill([baseline_history.generation; flipud(baseline_history.generation)], ...
     [baseline_history.best_obj; flipud(baseline_history.worst_obj)], ...
     COLOR_SCHEME.blue, 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
     'DisplayName', 'Solution range');

% Average solution
plot(baseline_history.generation, baseline_history.avg_obj, '--', ...
     'Color', COLOR_SCHEME.red, 'LineWidth', 3, 'DisplayName', 'Average');

% Best solution
plot(baseline_history.generation, baseline_history.best_obj, '-', ...
     'Color', COLOR_SCHEME.blue, 'LineWidth', 2.5, 'DisplayName', 'Best');

% Mark maximum improvement point
[max_improve, max_improve_gen] = max(diff(baseline_history.best_obj));
if max_improve < 0
    plot(baseline_history.generation(max_improve_gen), ...
         baseline_history.best_obj(max_improve_gen), ...
         'o', 'MarkerSize', 12, 'MarkerFaceColor', COLOR_SCHEME.red, ...
         'MarkerEdgeColor', 'w', 'LineWidth', 2, 'DisplayName', 'Max improvement');
end

% Mark convergence point (improvement rate < 0.1%)
improvement_rate = abs(diff(baseline_history.best_obj)) ./ baseline_history.best_obj(1:end-1) * 100;
convergence_gen = find(improvement_rate < 0.1, 1);
if ~isempty(convergence_gen)
    plot([convergence_gen convergence_gen], ylim, 'k--', 'LineWidth', 1, ...
         'HandleVisibility', 'off');
    text(convergence_gen, mean(ylim), sprintf(' Convergence\n (Gen %d)', convergence_gen), ...
         'FontSize', FONT_SIZE_LEGEND, 'Color', 'k');
end

xlabel('Generation', 'FontSize', FONT_SIZE_LABEL);
ylabel('Objective value', 'FontSize', FONT_SIZE_LABEL);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

% Performance info box
final_obj = baseline_history.best_obj(end);
initial_obj = baseline_history.best_obj(1);
total_improvement = (initial_obj - final_obj) / initial_obj * 100;
avg_improvement_per_gen = total_improvement / length(baseline_history.generation);

info_text = sprintf(['Initial: %.0f\n' ...
                     'Final: %.0f\n' ...
                     'Improvement: %.1f%%\n' ...
                     'Avg/Gen: %.2f%%\n' ...
                     'Converged: Gen %d'], ...
                    initial_obj, final_obj, total_improvement, ...
                    avg_improvement_per_gen, convergence_gen);
text(0.98, 0.98, info_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
     'FontSize', FONT_SIZE_LEGEND, 'BackgroundColor', 'w', ...
     'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 5);

fig2_path = fullfile(figures_dir, 'Fig2_Baseline_Convergence');
savefig(gcf, [fig2_path '.fig']);
fprintf('  Saved: Fig2_Baseline_Convergence.fig\n');

%% =========================================================================
% Figure 3: Baseline Route Map
% =========================================================================
fprintf('\nGenerating Figure 3: Baseline Routes...\n');
figure('Position', [150, 50, 1000, 750], 'Color', 'w', 'PaperPositionMode', 'auto');
hold on; box on;

baseline_colors = generate_academic_colors(baseline_NV);

% Draw routes with gradient effect
for i = 1:baseline_NV
    part_seq = baseline_VC{i};
    len = length(part_seq);
    route_color = baseline_colors(i,:);
    
    % Calculate total route length for gradient
    route_length = 0;
    for j = 1:len
        if j == 1
            route_length = route_length + norm(customer(part_seq(j),:) - vertexs(1,:));
        else
            route_length = route_length + norm(customer(part_seq(j),:) - customer(part_seq(j-1),:));
        end
    end
    route_length = route_length + norm(customer(part_seq(len),:) - vertexs(1,:));
    
    cum_length = 0;
    for j = 1:len
        if j == 1
            c1 = customer(part_seq(j),:);
            seg_length = norm(c1 - vertexs(1,:));
            alpha = cum_length / route_length;
            line_color = route_color * (1 - alpha * 0.5);
            plot([vertexs(1,1), c1(1)], [vertexs(1,2), c1(2)], ...
                 '-', 'Color', line_color, 'LineWidth', 2.5);
            
            % Direction arrow
            mid_x = (vertexs(1,1) + c1(1)) / 2;
            mid_y = (vertexs(1,2) + c1(2)) / 2;
            dx = c1(1) - vertexs(1,1);
            dy = c1(2) - vertexs(1,2);
            quiver(mid_x, mid_y, dx*0.1, dy*0.1, 0, 'Color', line_color, ...
                   'LineWidth', 1.5, 'MaxHeadSize', 1.5);
            
            cum_length = cum_length + seg_length;
        else
            c_pre = customer(part_seq(j-1),:);
            c_cur = customer(part_seq(j),:);
            seg_length = norm(c_cur - c_pre);
            alpha = cum_length / route_length;
            line_color = route_color * (1 - alpha * 0.5);
            plot([c_pre(1), c_cur(1)], [c_pre(2), c_cur(2)], ...
                 '-', 'Color', line_color, 'LineWidth', 2.5);
            
            % Direction arrow
            mid_x = (c_pre(1) + c_cur(1)) / 2;
            mid_y = (c_pre(2) + c_cur(2)) / 2;
            dx = c_cur(1) - c_pre(1);
            dy = c_cur(2) - c_pre(2);
            quiver(mid_x, mid_y, dx*0.1, dy*0.1, 0, 'Color', line_color, ...
                   'LineWidth', 1.5, 'MaxHeadSize', 1.5);
            
            cum_length = cum_length + seg_length;
        end
    end
    
    % Return line
    c_last = customer(part_seq(len),:);
    alpha = cum_length / route_length;
    line_color = route_color * (1 - alpha * 0.5);
    plot([c_last(1), vertexs(1,1)], [c_last(2), vertexs(1,2)], ...
         '--', 'Color', line_color, 'LineWidth', 2);
end

% Stations with bubble size = demand
bubble_sizes = 50 + abs(demands(:)) .* 15;
scatter(customer(:,1), customer(:,2), bubble_sizes, COLOR_SCHEME.blue, 'filled', ...
        'MarkerEdgeColor', 'b', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);

% Depot
plot(vertexs(1,1), vertexs(1,2), 'p', 'MarkerSize', 20, ...
     'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k', 'LineWidth', 2);

% Station ID labels
for i = 1:cusnum
    text(customer(i,1), customer(i,2), sprintf('%d', i), ...
         'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'EdgeColor', 'none');
end

xlabel('X coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Y coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);
axis equal tight;

% Info box
avg_load_rate = mean(calculate_load_rates(baseline_VC, demands, cap)) * 100;
info_text = sprintf(['Vehicles: %d\n' ...
                     'Total Distance: %.1f km\n' ...
                     'Total Cost: ¥%.0f\n' ...
                     'Avg Load Rate: %.1f%%\n' ...
                     'Violated Routes: %d'], ...
                    baseline_NV, baseline_TD, baseline_total_cost, ...
                    avg_load_rate, baseline_violate_num);
text(0.02, 0.98, info_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', FONT_SIZE_LEGEND, ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 5);

text(0.98, 0.02, 'Bubble size = Demand', 'Units', 'normalized', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
     'FontSize', FONT_SIZE_LEGEND, 'BackgroundColor', 'w', ...
     'EdgeColor', 'k', 'LineWidth', 1, 'Margin', 3);

fig3_path = fullfile(figures_dir, 'Fig3_Baseline_Routes');
savefig(gcf, [fig3_path '.fig']);
fprintf('  Saved: Fig3_Baseline_Routes.fig\n');

%% =========================================================================
% Step 2: Station Selection and User Dispatch Pairing
% =========================================================================
fprintf('\nStep 2: Station Analysis and User Dispatch Selection\n');
fprintf('-----------------------------------------------------------------\n');

dist_from_depot = dist(1, 2:end)';
remoteness_score = dist_from_depot / max(dist_from_depot);
demand_score = demands / max(demands);
low_demand_score = 1 - demand_score;
user_suitability_score = 0.7 * remoteness_score + 0.3 * low_demand_score;
suitable_mask = (demands <= 5);

% Calculate operator cost per station
operator_cost_per_station = zeros(cusnum, 1);
for i = 1:cusnum
    round_trip = 2 * dist_from_depot(i);
    travel_time = round_trip / avg_speed;
    operator_cost_per_station(i) = variable_cost_per_km * round_trip + ...
                                    driver_hourly_cost * travel_time + ...
                                    service_cost_per_station + ...
                                    (fixed_cost_per_vehicle * demands(i) / cap);
end

% Supply-demand analysis
rng(42);
supply = demands + randi([-5, 10], cusnum, 1);
demand_actual = demands + randi([-3, 8], cusnum, 1);
supply_demand_gap = supply - demand_actual;

% Calculate user dispatch incentive
user_incentive_per_station = zeros(cusnum, 1);
target_station_for_user = zeros(cusnum, 1);
bikes_to_dispatch = zeros(cusnum, 1);

for i = 1:cusnum
    best_target = 0;
    max_score = -inf;
    for j = 1:cusnum
        if i ~= j
            demand_score_val = max(0, -supply_demand_gap(j)) / max(abs(supply_demand_gap));
            distance_score = 1 / (1 + dist(i+1, j+1) / max(max(dist)));
            total_score = 0.7 * demand_score_val + 0.3 * distance_score;
            if total_score > max_score
                max_score = total_score;
                best_target = j;
            end
        end
    end
    target_station_for_user(i) = best_target;
    if best_target > 0
        actual_bikes = min([demands(i), 5]);
        bikes_to_dispatch(i) = actual_bikes;
        distance_to_target = dist(i+1, best_target+1);
        base_incentive = 5;
        distance_incentive = distance_to_target * 0.8;
        bikes_incentive = actual_bikes * 2.5;
        user_incentive_per_station(i) = base_incentive + distance_incentive + bikes_incentive;
        user_incentive_per_station(i) = max(8, min(50, user_incentive_per_station(i)));
    else
        user_incentive_per_station(i) = inf;
        bikes_to_dispatch(i) = 0;
    end
end

potential_saving = operator_cost_per_station - user_incentive_per_station;

% Select user dispatch stations
condition1 = suitable_mask;
condition2 = potential_saving > 0;
condition3 = user_incentive_per_station < inf;
valid_stations = find(condition1 & condition2 & condition3);
num_valid = length(valid_stations);

fprintf('  Valid candidate stations: %d\n', num_valid);

if num_valid == 0
    user_schedule_flag = zeros(cusnum, 1);
    num_user_stations = 0;
    total_user_incentive = 0;
    user_matching = [];
else
    [~, sorted_idx_in_valid] = sort(user_suitability_score(valid_stations), 'descend');
    sorted_valid_stations = valid_stations(sorted_idx_in_valid);
    
    % User station quantity control parameters
    target_user_ratio = 0.30;  % Target: 30%
    min_user_stations = 8;
    max_user_stations_ratio = floor(cusnum * target_user_ratio);
    max_user_stations = 15;  % Hard cap
    max_user_stations = min(max_user_stations, max_user_stations_ratio);
    min_operator_stations = 20;  % Ensure sufficient operator capacity
    
    num_user_stations = min([max_user_stations, num_valid, cusnum - min_operator_stations]);
    num_user_stations = max(min_user_stations, num_user_stations);
    num_user_stations = min(num_user_stations, num_valid);
    selected_stations = sorted_valid_stations(1:num_user_stations);
    user_schedule_flag = zeros(cusnum, 1);
    user_schedule_flag(selected_stations) = 1;
    total_user_incentive = sum(user_incentive_per_station(selected_stations));
    
    fprintf('  Selected user dispatch stations: %d (%.1f%%)\n', num_user_stations, num_user_stations/cusnum*100);
    
    user_matching = [];
    for i = 1:num_user_stations
        source_idx = selected_stations(i);
        target_idx = target_station_for_user(source_idx);
        if target_idx > 0
            match.source = source_idx;
            match.target = target_idx;
            match.bikes = bikes_to_dispatch(source_idx);
            match.distance = dist(source_idx+1, target_idx+1);
            match.incentive = user_incentive_per_station(source_idx);
            match.source_gap = supply_demand_gap(source_idx);
            match.target_gap = supply_demand_gap(target_idx);
            if isempty(user_matching)
                user_matching = match;
            else
                user_matching(end+1) = match;
            end
        end
    end
end

% =========================================================================
% Step 3: Hybrid Dispatch Optimization
% =========================================================================
fprintf('\nStep 3: Hybrid Dispatch Optimization\n');
fprintf('-----------------------------------------------------------------\n');

if num_user_stations == 0
    fprintf('  No user stations, using baseline solution\n');
    final_total_cost = baseline_total_cost;
    bestVC = baseline_VC;
    bestNV = baseline_NV;
    bestTD = baseline_TD;
    best_violate_num = baseline_violate_num;
    best_violate_cus = baseline_violate_cus;
    operator_total_cost = baseline_total_cost;
    hybrid_history = baseline_history;
    operator_customers = 1:cusnum;
    selected_stations = [];
    operator_vehicle_cost = baseline_vehicle_cost;
    operator_distance_cost = baseline_distance_cost;
    operator_time_cost = baseline_time_cost;
    operator_service_cost = baseline_service_cost;
    total_user_incentive = 0;
else
    operator_customers = find(~user_schedule_flag);
    operator_cusnum = length(operator_customers);
    fprintf('  Operator dispatch stations: %d\n', operator_cusnum);
    
    N = operator_cusnum + v_num - 1;
    depot_and_operator = [1; operator_customers+1];
    sub_dist = dist(depot_and_operator, depot_and_operator);
    init_vc = init(operator_cusnum, a(operator_customers), demands(operator_customers), cap);
    Chrom = InitPopCW(NIND, N, operator_cusnum, init_vc);
    
    % Record hybrid iteration history
    MAXGEN_HYBRID = 100;
    hybrid_history = struct();
    hybrid_history.generation = zeros(MAXGEN_HYBRID, 1);
    hybrid_history.best_obj = zeros(MAXGEN_HYBRID, 1);
    hybrid_history.avg_obj = zeros(MAXGEN_HYBRID, 1);
    hybrid_history.worst_obj = zeros(MAXGEN_HYBRID, 1);
    
    fprintf('  Hybrid optimization in progress...');
    gen = 1;
    while gen <= MAXGEN_HYBRID
        ObjV = calObj_hybrid(Chrom, operator_cusnum, cap, ...
            demands(operator_customers), a(operator_customers), ...
            b(operator_customers), L, s(operator_customers), ...
            sub_dist, alpha, beta, 0, 1, operator_customers);
        
        hybrid_history.generation(gen) = gen;
        hybrid_history.best_obj(gen) = min(ObjV);
        hybrid_history.avg_obj(gen) = mean(ObjV);
        hybrid_history.worst_obj(gen) = max(ObjV);
        
        FitnV = Fitness(ObjV);
        SelCh = Select(Chrom, FitnV, GGAP);
        SelCh = Recombin(SelCh, Pc);
        SelCh = Mutate(SelCh, Pm);
        SelCh = LocalSearch_hybrid(SelCh, operator_cusnum, cap, ...
            demands(operator_customers), a(operator_customers), ...
            b(operator_customers), L, s(operator_customers), ...
            sub_dist, alpha, beta, operator_customers);
        Chrom = Reins(Chrom, SelCh, ObjV);
        Chrom = deal_Repeat(Chrom);
        
        if mod(gen, 20) == 0
            fprintf('.');
        end
        gen = gen + 1;
    end
    fprintf('Complete\n');
    
    ObjV = calObj_hybrid(Chrom, operator_cusnum, cap, ...
        demands(operator_customers), a(operator_customers), ...
        b(operator_customers), L, s(operator_customers), ...
        sub_dist, alpha, beta, 0, 1, operator_customers);
    [~, minInd] = min(ObjV);
    [bestVC, bestNV, bestTD, best_violate_num, best_violate_cus] = decode_hybrid(...
        Chrom(minInd(1),:), operator_cusnum, cap, demands(operator_customers), ...
        a(operator_customers), b(operator_customers), L, ...
        s(operator_customers), sub_dist, operator_customers);
    
    operator_time_hours = bestTD / avg_speed;
    operator_vehicle_cost = bestNV * fixed_cost_per_vehicle;
    operator_distance_cost = bestTD * variable_cost_per_km;
    operator_time_cost = operator_time_hours * driver_hourly_cost;
    operator_service_cost = operator_cusnum * service_cost_per_station;
    operator_total_cost = operator_vehicle_cost + operator_distance_cost + ...
                          operator_time_cost + operator_service_cost;
    final_total_cost = operator_total_cost + total_user_incentive;
end

fprintf('\nHybrid Solution Results:\n');
fprintf('  Vehicles: %d\n', bestNV);
fprintf('  Total Distance: %.2f km\n', bestTD);
fprintf('  Violated Routes: %d\n', best_violate_num);
fprintf('  Violated Stations: %d\n', best_violate_cus);
fprintf('  Operator Cost: ¥%.2f\n', operator_total_cost);
fprintf('  User Incentive: ¥%.2f\n', total_user_incentive);
fprintf('  Total Cost: ¥%.2f\n', final_total_cost);

cost_diff = baseline_total_cost - final_total_cost;
fprintf('\nCost Savings: ¥%.2f (%.2f%%)\n', cost_diff, cost_diff/baseline_total_cost*100);

%% =========================================================================
% Figure 4: Hybrid Algorithm Convergence Curve
% =========================================================================
fprintf('\nGenerating Figure 4: Hybrid Convergence...\n');
figure('Position', [200, 50, 1000, 750], 'Color', 'w', 'PaperPositionMode', 'auto');
hold on; box on;

% Fill hybrid solution range
fill([hybrid_history.generation; flipud(hybrid_history.generation)], ...
     [hybrid_history.best_obj; flipud(hybrid_history.worst_obj)], ...
     COLOR_SCHEME.green, 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
     'DisplayName', 'Hybrid solution range');

% Baseline best as reference
max_gen = min(length(baseline_history.generation), length(hybrid_history.generation));
plot(baseline_history.generation(1:max_gen), baseline_history.best_obj(1:max_gen), '--', ...
     'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', 'Baseline best');

% Hybrid average
plot(hybrid_history.generation, hybrid_history.avg_obj, ':', ...
     'Color', COLOR_SCHEME.orange, 'LineWidth', 1.5, 'DisplayName', 'Hybrid average');

% Hybrid best
plot(hybrid_history.generation, hybrid_history.best_obj, '-', ...
     'Color', COLOR_SCHEME.green, 'LineWidth', 2.5, 'DisplayName', 'Hybrid best');

% Mark improvement percentage
if num_user_stations > 0
    baseline_final = baseline_history.best_obj(max_gen);
    hybrid_final = hybrid_history.best_obj(max_gen);
    improvement_pct = (baseline_final - hybrid_final) / baseline_final * 100;
    
    y_pos = (baseline_final + hybrid_final) / 2;
    plot([max_gen, max_gen], [hybrid_final, baseline_final], 'r-', 'LineWidth', 2);
    text(max_gen + 2, y_pos, sprintf('%.1f%% improvement', improvement_pct), ...
         'FontSize', FONT_SIZE_LEGEND, 'Color', 'r', 'FontWeight', 'bold', ...
         'BackgroundColor', 'w', 'EdgeColor', 'r', 'Margin', 2);
end

% Mark key iteration points
[max_improve, max_improve_gen] = max(diff(hybrid_history.best_obj));
if max_improve < 0
    plot(hybrid_history.generation(max_improve_gen), ...
         hybrid_history.best_obj(max_improve_gen), ...
         'o', 'MarkerSize', 12, 'MarkerFaceColor', COLOR_SCHEME.red, ...
         'MarkerEdgeColor', 'w', 'LineWidth', 2);
end

% Mark convergence
improvement_rate = abs(diff(hybrid_history.best_obj)) ./ hybrid_history.best_obj(1:end-1) * 100;
convergence_gen = find(improvement_rate < 0.1, 1);
if ~isempty(convergence_gen)
    plot([convergence_gen convergence_gen], ylim, 'k--', 'LineWidth', 1, ...
         'HandleVisibility', 'off');
    text(convergence_gen, mean(ylim), sprintf(' Convergence\n (Gen %d)', convergence_gen), ...
         'FontSize', FONT_SIZE_LEGEND-1, 'Color', 'k');
end

xlabel('Generation', 'FontSize', FONT_SIZE_LABEL);
ylabel('Objective value', 'FontSize', FONT_SIZE_LABEL);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

% Performance comparison info box
if num_user_stations > 0
    hybrid_initial = hybrid_history.best_obj(1);
    hybrid_final_val = hybrid_history.best_obj(end);
    hybrid_improvement = (hybrid_initial - hybrid_final_val) / hybrid_initial * 100;
    
    info_text = sprintf(['Hybrid Algorithm:\n' ...
                         'Initial: %.0f\n' ...
                         'Final: %.0f\n' ...
                         'Improvement: %.1f%%\n' ...
                         'vs Baseline: %.1f%% better'], ...
                        hybrid_initial, hybrid_final_val, hybrid_improvement, improvement_pct);
else
    info_text = sprintf(['Initial: %.0f\n' ...
                         'Final: %.0f\n' ...
                         'Improvement: %.1f%%'], ...
                        hybrid_history.best_obj(1), hybrid_history.best_obj(end), ...
                        (hybrid_history.best_obj(1) - hybrid_history.best_obj(end)) / ...
                        hybrid_history.best_obj(1) * 100);
end
text(0.98, 0.98, info_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
     'FontSize', FONT_SIZE_LEGEND, 'BackgroundColor', 'w', ...
     'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 5);

fig4_path = fullfile(figures_dir, 'Fig4_Hybrid_Convergence');
savefig(gcf, [fig4_path '.fig']);
fprintf('  Saved: Fig4_Hybrid_Convergence.fig\n');

%% =========================================================================
% Figure 5: Hybrid Route Map
% =========================================================================
fprintf('\nGenerating Figure 5: Hybrid Routes...\n');
figure('Position', [250, 50, 1000, 750], 'Color', 'w', 'PaperPositionMode', 'auto');
hold on; box on;

if num_user_stations > 0
    hybrid_colors = generate_academic_colors(bestNV);
    
    % Draw routes with gradient and direction arrows
    for i = 1:bestNV
        part_seq = bestVC{i};
        len = length(part_seq);
        actual_idx = operator_customers(part_seq);
        route_color = hybrid_colors(i,:);
        
        % Calculate route length
        route_length = 0;
        for j = 1:len
            if j == 1
                route_length = route_length + norm(customer(actual_idx(j),:) - vertexs(1,:));
            else
                route_length = route_length + norm(customer(actual_idx(j),:) - customer(actual_idx(j-1),:));
            end
        end
        route_length = route_length + norm(customer(actual_idx(len),:) - vertexs(1,:));
        
        cum_length = 0;
        for j = 1:len
            if j == 1
                c1 = customer(actual_idx(j),:);
                seg_length = norm(c1 - vertexs(1,:));
                alpha = cum_length / route_length;
                line_color = route_color * (1 - alpha * 0.5);
                plot([vertexs(1,1), c1(1)], [vertexs(1,2), c1(2)], ...
                     '-', 'Color', line_color, 'LineWidth', 2.5);
                     
                % Direction arrow
                mid_x = (vertexs(1,1) + c1(1)) / 2;
                mid_y = (vertexs(1,2) + c1(2)) / 2;
                dx = c1(1) - vertexs(1,1);
                dy = c1(2) - vertexs(1,2);
                quiver(mid_x, mid_y, dx*0.1, dy*0.1, 0, 'Color', line_color, ...
                       'LineWidth', 1.5, 'MaxHeadSize', 1.5);
                
                cum_length = cum_length + seg_length;
            else
                c_pre = customer(actual_idx(j-1),:);
                c_cur = customer(actual_idx(j),:);
                seg_length = norm(c_cur - c_pre);
                alpha = cum_length / route_length;
                line_color = route_color * (1 - alpha * 0.5);
                plot([c_pre(1), c_cur(1)], [c_pre(2), c_cur(2)], ...
                     '-', 'Color', line_color, 'LineWidth', 2.5);
                     
                % Direction arrow
                mid_x = (c_pre(1) + c_cur(1)) / 2;
                mid_y = (c_pre(2) + c_cur(2)) / 2;
                dx = c_cur(1) - c_pre(1);
                dy = c_cur(2) - c_pre(2);
                quiver(mid_x, mid_y, dx*0.1, dy*0.1, 0, 'Color', line_color, ...
                       'LineWidth', 1.5, 'MaxHeadSize', 1.5);
                
                cum_length = cum_length + seg_length;
            end
        end
        
        % Return line
        c_last = customer(actual_idx(len),:);
        alpha = cum_length / route_length;
        line_color = route_color * (1 - alpha * 0.5);
        plot([c_last(1), vertexs(1,1)], [c_last(2), vertexs(1,2)], ...
             '--', 'Color', line_color, 'LineWidth', 2);
    end
    
    % Operator stations (bubble size = demand)
    operator_station_coords = customer(operator_customers,:);
    operator_demands = demands(operator_customers);
    bubble_sizes_op = 50 + abs(operator_demands(:)) .* 10;
    scatter(operator_station_coords(:,1), operator_station_coords(:,2), ...
            bubble_sizes_op, COLOR_SCHEME.blue, 'o', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
    
    % Operator station labels
    for i = 1:length(operator_customers)
        text(operator_station_coords(i,1), operator_station_coords(i,2), ...
             sprintf('%d', operator_customers(i)), ...
             'FontSize', 7, 'Color', 'k', 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'EdgeColor', 'none');
    end
    
    % User stations (large triangle markers)
    user_station_idx = find(user_schedule_flag);
    user_station_coords = customer(user_station_idx,:);
    scatter(user_station_coords(:,1), user_station_coords(:,2), 250, ...
            COLOR_SCHEME.green, '^', 'filled', ...
            'LineWidth', 2.5);
    
    % User station labels
    for i = 1:length(user_station_idx)
        text(user_station_coords(i,1), user_station_coords(i,2), ...
             sprintf('%d', user_station_idx(i)), ...
             'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'EdgeColor', 'none');
    end
    
    % Depot
    plot(vertexs(1,1), vertexs(1,2), 'p', 'MarkerSize', 20, ...
         'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Legend
    legend_handles = [
        scatter(NaN, NaN, 100, COLOR_SCHEME.blue, 'o', 'filled', 'MarkerEdgeColor', 'w'), ...
        scatter(NaN, NaN, 250, COLOR_SCHEME.green, '^', 'filled', 'LineWidth', 2), ...
        plot(NaN, NaN, 'p', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k')
    ];
    legend(legend_handles, {'Operator stations', 'User stations (highlighted)', 'Depot'}, ...
           'Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
else
    % Use baseline plot when no user stations
    hybrid_colors = generate_academic_colors(bestNV);
    for i = 1:bestNV
        part_seq = bestVC{i};
        len = length(part_seq);
        route_color = hybrid_colors(i,:);
        
        for j = 1:len
            if j == 1
                c1 = customer(part_seq(j),:);
                plot([vertexs(1,1), c1(1)], [vertexs(1,2), c1(2)], ...
                     '-', 'Color', route_color, 'LineWidth', LINE_WIDTH_MAIN);
            else
                c_pre = customer(part_seq(j-1),:);
                c_cur = customer(part_seq(j),:);
                plot([c_pre(1), c_cur(1)], [c_pre(2), c_cur(2)], ...
                     '-', 'Color', route_color, 'LineWidth', LINE_WIDTH_MAIN);
            end
        end
        c_last = customer(part_seq(len),:);
        plot([c_last(1), vertexs(1,1)], [c_last(2), vertexs(1,2)], ...
             '-', 'Color', route_color, 'LineWidth', LINE_WIDTH_MAIN);
    end
    
    plot(customer(:,1), customer(:,2), 'o', 'MarkerSize', MARKER_SIZE, ...
         'MarkerFaceColor', COLOR_SCHEME.blue, 'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
    plot(vertexs(1,1), vertexs(1,2), 's', 'MarkerSize', 12, ...
         'MarkerFaceColor', COLOR_SCHEME.depot, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

xlabel('X coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Y coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);
axis equal tight;

% Detailed info box
if num_user_stations > 0
    cost_saving = baseline_total_cost - final_total_cost;
    cost_saving_pct = cost_saving / baseline_total_cost * 100;
    avg_load_hybrid = mean(calculate_load_rates(bestVC, demands(operator_customers), cap)) * 100;
    
    info_text = sprintf(['Vehicles: %d (-%d)\n' ...
                         'Distance: %.1f km (%.1f%%)\n' ...
                         'User Stations: %d\n' ...
                         'Total Cost: ¥%.0f\n' ...
                         'Cost Saving: ¥%.0f (%.1f%%)\n' ...
                         'Avg Load: %.1f%%'], ...
                        bestNV, baseline_NV - bestNV, bestTD, ...
                        (baseline_TD - bestTD)/baseline_TD*100, ...
                        num_user_stations, final_total_cost, ...
                        cost_saving, cost_saving_pct, avg_load_hybrid);
else
    info_text = sprintf('Vehicles: %d\nDistance: %.1f km\nCost: ¥%.0f', ...
                        bestNV, bestTD, final_total_cost);
end

text(0.02, 0.98, info_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', FONT_SIZE_LEGEND, ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 5);

if num_user_stations > 0
    text(0.98, 0.02, 'Bubble size = Operator demand', 'Units', 'normalized', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
         'FontSize', FONT_SIZE_LEGEND, 'BackgroundColor', 'w', ...
         'EdgeColor', 'k', 'LineWidth', 1, 'Margin', 3);
end

fig5_path = fullfile(figures_dir, 'Fig5_Hybrid_Routes');
savefig(gcf, [fig5_path '.fig']);
fprintf('  Saved: Fig5_Hybrid_Routes.fig\n');

%% =========================================================================
% Figure 6: Incentive Coefficient Comparison
% =========================================================================
fprintf('\nGenerating Figure 6: Incentive Comparison...\n');

if num_user_stations > 0
    % Calculate costs under different incentive coefficients
    incentive_coefficients = [0.7, 0.85, 1.0];
    num_coeff = length(incentive_coefficients);
    incentive_costs_over_time = zeros(num_user_stations, num_coeff);
    
    for c = 1:num_coeff
        coeff = incentive_coefficients(c);
        for i = 1:num_user_stations
            station_idx = selected_stations(i);
            incentive_costs_over_time(i, c) = user_incentive_per_station(station_idx) * coeff;
        end
    end
    
    figure('Position', [300, 50, 1000, 750], 'Color', 'w', 'PaperPositionMode', 'auto');
    hold on; box on;
    
    colors_coeff = {COLOR_SCHEME.blue, COLOR_SCHEME.orange, COLOR_SCHEME.red};
    line_styles = {'-', '--', ':'};
    
    % Cumulative cost curves
    for c = 1:num_coeff
        cumulative_cost = cumsum(incentive_costs_over_time(:, c));
        plot(1:num_user_stations, cumulative_cost, line_styles{c}, ...
             'Color', colors_coeff{c}, 'LineWidth', 3, ...
             'DisplayName', sprintf('λ = %.2f', incentive_coefficients(c)));
    end
    
    % Marginal cost analysis
    yyaxis right;
    marginal_cost_085 = [incentive_costs_over_time(1, 2); diff(cumsum(incentive_costs_over_time(:, 2)))];
    plot(1:num_user_stations, marginal_cost_085, 'o-', ...
         'Color', [1 0.1 0.1], 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'MarkerFaceColor', [0.7 0.7 0.7], 'DisplayName', 'Marginal cost (λ=0.85)');
    ylabel('Marginal cost per station (CNY)', 'FontSize', FONT_SIZE_LABEL);
    set(gca, 'YColor', [0.5 0.5 0.5]);
    
    yyaxis left;
    ylabel('Cumulative incentive cost (CNY)', 'FontSize', FONT_SIZE_LABEL);
    set(gca, 'YColor', 'k');
    
    % Highlight recommended coefficient
    final_costs = cumsum(incentive_costs_over_time(:, :));
    [~, optimal_idx] = min(final_costs(end, :));
    optimal_coeff = incentive_coefficients(optimal_idx);
    
    if optimal_coeff == 0.85
        y_lim = ylim;
        fill([0.5 num_user_stations+0.5 num_user_stations+0.5 0.5], ...
             [y_lim(1) y_lim(1) y_lim(2) y_lim(2)], ...
             COLOR_SCHEME.orange, 'FaceAlpha', 0.05, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
    end
    
    xlabel('Number of user-dispatched stations', 'FontSize', FONT_SIZE_LABEL);
    legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
    grid on;
    set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);
        
    fig6_path = fullfile(figures_dir, 'Fig6_Incentive_Comparison');
    savefig(gcf, [fig6_path '.fig']);
    fprintf('  Saved: Fig6_Incentive_Comparison.fig\n');
else
    fprintf('  Skipped Figure 6 (no user stations)\n');
end

%% =========================================================================
% Figure 7: Cumulative Station Weight
% =========================================================================
fprintf('\nGenerating Figure 7: Cumulative Weight...\n');

if num_user_stations > 0
    % Calculate node weight (based on suitability score)
    selected_scores = user_suitability_score(selected_stations);
    [sorted_scores, sort_idx] = sort(selected_scores, 'descend');
    cumulative_weight = cumsum(sorted_scores) / sum(sorted_scores);
    
    figure('Position', [350, 50, 1000, 750], 'Color', 'w', 'PaperPositionMode', 'auto');
    hold on; box on;
    
    % Cumulative curve with area fill
    area(1:num_user_stations, cumulative_weight * 100, ...
         'FaceColor', COLOR_SCHEME.blue, 'FaceAlpha', 0.2, ...
         'EdgeColor', COLOR_SCHEME.blue, 'LineWidth', 2.5);
    
    % Marker points
    plot(1:num_user_stations, cumulative_weight * 100, ...
         'o', 'MarkerSize', 7, 'MarkerFaceColor', COLOR_SCHEME.blue, ...
         'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    
    % Multiple threshold reference lines
    thresholds = [70, 80, 90];
    threshold_colors = {[0.7 0.7 0.7], [0.5 0.5 0.5], [0.3 0.3 0.3]};
    threshold_stations = zeros(length(thresholds), 1);
    
    for t = 1:length(thresholds)
        thresh = thresholds(t);
        thresh_idx = find(cumulative_weight * 100 >= thresh, 1);
        if ~isempty(thresh_idx)
            threshold_stations(t) = thresh_idx;
            
            % Horizontal line
            plot([1, num_user_stations], [thresh, thresh], '--', ...
                 'Color', threshold_colors{t}, 'LineWidth', 1.5);
            
            % Vertical line
            plot([thresh_idx, thresh_idx], [0, thresh], ':', ...
                 'Color', threshold_colors{t}, 'LineWidth', 1.5);
            
            % Marker point
            plot(thresh_idx, thresh, 's', 'MarkerSize', 10, ...
                 'MarkerFaceColor', COLOR_SCHEME.red, ...
                 'MarkerEdgeColor', 'w', 'LineWidth', 2);
            
            % Label text
            text(num_user_stations * 0.85, thresh + 2, ...
                 sprintf('%d%% (%d stations)', thresh, thresh_idx), ...
                 'FontSize', FONT_SIZE_LEGEND, 'Color', threshold_colors{t}, ...
                 'FontWeight', 'bold', 'BackgroundColor', 'w');
        end
    end
    
    % Mark Pareto frontier (most efficient first few stations)
    pareto_count = min(5, num_user_stations);
    for i = 1:pareto_count
        marginal_contrib = 0;
        if i == 1
            marginal_contrib = cumulative_weight(i) * 100;
        else
            marginal_contrib = (cumulative_weight(i) - cumulative_weight(i-1)) * 100;
        end
        
        if marginal_contrib > 10  % Only mark significant contributions
            text(i, cumulative_weight(i) * 100 - 5, ...
                 sprintf('%.1f%%', marginal_contrib), ...
                 'FontSize', FONT_SIZE_LEGEND-1, 'Color', COLOR_SCHEME.red, ...
                 'HorizontalAlignment', 'center');
        end
    end
    
    xlabel('Number of user-dispatched stations', 'FontSize', FONT_SIZE_LABEL);
    ylabel('Cumulative weight (%)', 'FontSize', FONT_SIZE_LABEL);
    grid on;
    set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);
    ylim([0, 105]);
    xlim([0.5, num_user_stations + 0.5]);
    
    % Strategy info box
    info_text = sprintf(['Selection Strategy:\n' ...
                         '70%% coverage: %d stations\n' ...
                         '80%% coverage: %d stations ✓\n' ...
                         '90%% coverage: %d stations\n' ...
                         'Efficiency: %.1f%%/station'], ...
                        threshold_stations(1), threshold_stations(2), ...
                        threshold_stations(3), ...
                        80 / threshold_stations(2));
    text(0.02, 0.98, info_text, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'FontSize', FONT_SIZE_LEGEND, ...
         'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Margin', 5);
    
    fig7_path = fullfile(figures_dir, 'Fig7_Cumulative_Weight');
    savefig(gcf, [fig7_path '.fig']);
    fprintf('  Saved: Fig7_Cumulative_Weight.fig\n');
else
    fprintf('  Skipped Figure 7 (no user stations)\n');
end

%% =========================================================================
% Figure 8: User Dispatch Network Analysis
% =========================================================================
if num_user_stations > 0 && ~isempty(user_matching)
    fprintf('\nGenerating Figure 8: User Dispatch Network...\n');
    figure('Position', [400, 50, 800, 650], 'Color', 'w', 'PaperPositionMode', 'auto');
    
    hold on; box on;
    
    % Voronoi background (faded)
    [vx, vy] = voronoi(customer(:,1), customer(:,2));
    plot(vx, vy, '-', 'Color', [0.95 0.95 0.95], 'LineWidth', 0.5);
    
    % All operator stations (small dots)
    operator_coords = customer(operator_customers, :);
    scatter(operator_coords(:,1), operator_coords(:,2), 100, ...
            'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.5);
    
    % Operator station labels
    for i = 1:length(operator_customers)
        text(operator_station_coords(i,1), operator_station_coords(i,2), ...
             sprintf('%d', operator_customers(i)), ...
             'FontSize', 7, 'Color', 'k', 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'EdgeColor', 'none');
    end

    % Matching lines (with gradient and arrows)
    for i = 1:length(user_matching)
        source_coords = customer(user_matching(i).source, :);
        target_coords = customer(user_matching(i).target, :);
        
        % Line width and color based on incentive
        incentive_normalized = (user_matching(i).incentive - 8) / (50 - 8);
        line_width = 1 + incentive_normalized * 2.5;
        line_color = [0.2 + incentive_normalized*0.5, 0.6 - incentive_normalized*0.4, 0.2];
        
        % Curved connection line
        t = linspace(0, 1, 20);
        mid_offset_x = (target_coords(1) - source_coords(1)) * 0.1 * sin(pi);
        mid_offset_y = (target_coords(2) - source_coords(2)) * 0.1 * sin(pi);
        
        curve_x = source_coords(1) + t.*(target_coords(1) - source_coords(1)) + ...
                  mid_offset_x * sin(pi*t);
        curve_y = source_coords(2) + t.*(target_coords(2) - source_coords(2)) + ...
                  mid_offset_y * sin(pi*t);
        
        plot(curve_x, curve_y, '-', 'Color', line_color, ...
             'LineWidth', line_width, 'LineStyle', '--');
        
        % Direction arrow
        arrow_idx = 15;
        dx = curve_x(arrow_idx+1) - curve_x(arrow_idx);
        dy = curve_y(arrow_idx+1) - curve_y(arrow_idx);
        quiver(curve_x(arrow_idx), curve_y(arrow_idx), dx*2, dy*2, 0, ...
               'Color', line_color, 'LineWidth', line_width, 'MaxHeadSize', 2);
    end
    
    % User stations (large triangle markers)
    user_coords = customer(selected_stations, :);
    scatter(user_coords(:,1), user_coords(:,2), 300, ...
            COLOR_SCHEME.green, '^', 'filled', ...
            'LineWidth', 2.5);
    
    % User station labels
    for i = 1:length(selected_stations)
        text(user_coords(i,1), user_coords(i,2), sprintf('%d', selected_stations(i)), ...
             'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'EdgeColor', 'none');
    end
    
    % Depot
    plot(vertexs(1,1), vertexs(1,2), 'p', 'MarkerSize', 18, ...
         'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    xlabel('X coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
    ylabel('Y coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
    grid on;
    set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);
    axis equal tight;
    xlim([115,122]);
    ylim([21,27]);

    % Legend
    legend_items = [
        scatter(NaN, NaN, 300, COLOR_SCHEME.green, '^', 'filled', 'LineWidth', 2), ...
        plot(NaN, NaN, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.7 0.7 0.7]), ...
        plot(NaN, NaN, '--', 'Color', [0.4 0.5 0.3], 'LineWidth', 2)
    ];
    legend(legend_items, {'User stations', 'Operator stations', 'Dispatch flow'}, ...
           'Location', 'best', 'FontSize', FONT_SIZE_LEGEND-1, 'Box', 'off');
end

%% =========================================================================
% Figure 9: Station Clustering and Demand Heatmap
% =========================================================================
figure('Position', [100, 100, 600, 400], 'Color', 'w', 'PaperPositionMode', 'auto');

hold on; box on;

% K-means clustering
rng(42);
num_clusters = 5;
[cluster_idx, cluster_centers] = kmeans(customer, num_clusters, 'Replicates', 10);

% Plot clustering results
cluster_colors = lines(num_clusters);
for i = 1:num_clusters
    cluster_points = customer(cluster_idx == i, :);
    scatter(cluster_points(:,1), cluster_points(:,2), 150, cluster_colors(i,:), ...
            'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, ...
            'DisplayName', sprintf('Cluster %d', i));
end

% Operator station labels
for i = 1:length(operator_customers)
    text(operator_station_coords(i,1), operator_station_coords(i,2), ...
        sprintf('%d', operator_customers(i)), ...
        'FontSize', 7, 'Color', 'k', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none');
end

% User station labels
for i = 1:length(selected_stations)
    text(user_coords(i,1), user_coords(i,2), sprintf('%d', selected_stations(i)), ...
        'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none');
end

% Cluster centers
scatter(cluster_centers(:,1), cluster_centers(:,2), 200, 'k', 'p', ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 2, ...
        'DisplayName', 'Centers');

% Depot
plot(vertexs(1,1), vertexs(1,2), 's', 'MarkerSize', 15, ...
     'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', ...
     'LineWidth', 2, 'DisplayName', 'Depot');
xlim([115,122]);
ylim([21,27]);
xlabel('X coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Y coordinate (km)', 'FontSize', FONT_SIZE_LABEL);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-1, 'Box', 'off', 'NumColumns', 2);
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);
axis equal tight;

%% Cluster demand distribution boxplot
figure('Position', [100, 100, 600, 400], 'Color', 'w', 'PaperPositionMode', 'auto');

hold on; box on;

cluster_demands = cell(num_clusters, 1);
for i = 1:num_clusters
    cluster_demands{i} = demands(cluster_idx == i);
end

% Boxplot
positions = 1:num_clusters;
for i = 1:num_clusters
    boxplot_data = cluster_demands{i};
    bp = boxplot(boxplot_data, 'Positions', i, 'Width', 0.6, ...
                 'Colors', cluster_colors(i,:), 'Symbol', '');
    
    % Jittered scatter points
    x_jitter = i + 0.15*randn(size(boxplot_data));
    scatter(x_jitter, boxplot_data, 20, cluster_colors(i,:), ...
            'filled', 'MarkerFaceAlpha', 0.4);
end

xlabel('Cluster ID', 'FontSize', FONT_SIZE_LABEL);
ylabel('Demand (bikes)', 'FontSize', FONT_SIZE_LABEL);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND-1, 'Box', 'off', 'NumColumns', 2);

grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);
xlim([0.5, num_clusters+0.5]);

%% =========================================================================
% Figure 10: 3D Route Efficiency and Load Analysis
% =========================================================================
fprintf('\nGenerating Figure 10: 3D Efficiency Analysis...\n');

% Calculate load rates
baseline_loads = calculate_load_rates(baseline_VC, demands, cap);
if num_user_stations > 0
    hybrid_loads = calculate_load_rates(bestVC, demands(operator_customers), cap);
else
    hybrid_loads = baseline_loads;
end

figure('Position', [150, 100, 1200, 900], 'Color', 'w', 'PaperPositionMode', 'auto');

% Subplot 1: 3D scatter - route length vs stations vs load rate
subplot(2, 2, 1);
hold on; grid on; box on;

% Baseline data
baseline_route_lengths_local = zeros(baseline_NV, 1);
baseline_route_stations = zeros(baseline_NV, 1);
for i = 1:baseline_NV
    baseline_route_lengths_local(i) = part_length(baseline_VC{i}, dist);
    baseline_route_stations(i) = length(baseline_VC{i});
end

scatter3(baseline_route_lengths_local, baseline_route_stations, baseline_loads*100, ...
         100, COLOR_SCHEME.blue, 'filled', 'MarkerEdgeColor', 'w', ...
         'LineWidth', 1.2, 'DisplayName', 'Baseline');

if num_user_stations > 0
    % Hybrid data
    hybrid_route_lengths_local = zeros(bestNV, 1);
    hybrid_route_stations = zeros(bestNV, 1);
    for i = 1:bestNV
        global_route = operator_customers(bestVC{i});
        hybrid_route_lengths_local(i) = part_length(global_route, dist);
        hybrid_route_stations(i) = length(bestVC{i});
    end
    
    scatter3(hybrid_route_lengths_local, hybrid_route_stations, hybrid_loads*100, ...
             100, COLOR_SCHEME.green, '^', 'filled', 'MarkerEdgeColor', 'w', ...
             'LineWidth', 1.2, 'DisplayName', 'Hybrid');
end

xlabel('Route length (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Number of stations', 'FontSize', FONT_SIZE_LABEL);
zlabel('Load rate (%)', 'FontSize', FONT_SIZE_LABEL);
title('(a) 3D route efficiency space', 'FontSize', FONT_SIZE_SUBPLOT);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND);
view(45, 25);
set(gca, 'GridAlpha', 0.3);

% Subplot 2: Load rate-distance relationship (with regression)
subplot(2, 2, 2);
hold on; box on;

% Baseline
scatter(baseline_route_lengths_local, baseline_loads*100, 80, ...
        COLOR_SCHEME.blue, 'filled', 'MarkerEdgeColor', 'w', ...
        'LineWidth', 1, 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Baseline routes');

% Fit curve
p_baseline = polyfit(baseline_route_lengths_local, baseline_loads*100, 2);
x_fit = linspace(min(baseline_route_lengths_local), max(baseline_route_lengths_local), 100);
y_fit_baseline = polyval(p_baseline, x_fit);
plot(x_fit, y_fit_baseline, '--', 'Color', COLOR_SCHEME.blue, ...
     'LineWidth', 2, 'DisplayName', 'Baseline trend');

if num_user_stations > 0
    scatter(hybrid_route_lengths_local, hybrid_loads*100, 80, ...
            COLOR_SCHEME.green, '^', 'filled', 'MarkerEdgeColor', 'w', ...
            'LineWidth', 1, 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Hybrid routes');
    
    p_hybrid = polyfit(hybrid_route_lengths_local, hybrid_loads*100, 2);
    y_fit_hybrid = polyval(p_hybrid, x_fit);
    plot(x_fit, y_fit_hybrid, '--', 'Color', COLOR_SCHEME.green, ...
         'LineWidth', 2, 'DisplayName', 'Hybrid trend');
end

xlabel('Route length (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Load rate (%)', 'FontSize', FONT_SIZE_LABEL);
title('(b) Load rate vs route length', 'FontSize', FONT_SIZE_SUBPLOT);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

% Subplot 3: Load rate distribution (kernel density estimation)
subplot(2, 2, 3);
hold on; box on;

% KDE
[f_baseline, xi_baseline] = ksdensity(baseline_loads*100, 'BoundaryCorrection', 'reflection');
area(xi_baseline, f_baseline, 'FaceColor', COLOR_SCHEME.blue, ...
     'FaceAlpha', 0.3, 'EdgeColor', COLOR_SCHEME.blue, 'LineWidth', 2, ...
     'DisplayName', 'Baseline');

if num_user_stations > 0
    [f_hybrid, xi_hybrid] = ksdensity(hybrid_loads*100, 'BoundaryCorrection', 'reflection');
    area(xi_hybrid, f_hybrid, 'FaceColor', COLOR_SCHEME.green, ...
         'FaceAlpha', 0.3, 'EdgeColor', COLOR_SCHEME.green, 'LineWidth', 2, ...
         'DisplayName', 'Hybrid');
end

% Mean lines
mean_baseline = mean(baseline_loads*100);
plot([mean_baseline mean_baseline], ylim, '--', 'Color', COLOR_SCHEME.blue, ...
     'LineWidth', 1.5, 'HandleVisibility', 'off');
text(mean_baseline, max(ylim)*0.9, sprintf('μ=%.1f%%', mean_baseline), ...
     'FontSize', FONT_SIZE_LEGEND, 'Color', COLOR_SCHEME.blue);

if num_user_stations > 0
    mean_hybrid = mean(hybrid_loads*100);
    plot([mean_hybrid mean_hybrid], ylim, '--', 'Color', COLOR_SCHEME.green, ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(mean_hybrid, max(ylim)*0.8, sprintf('μ=%.1f%%', mean_hybrid), ...
         'FontSize', FONT_SIZE_LEGEND, 'Color', COLOR_SCHEME.green);
end

xlabel('Load rate (%)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Probability density', 'FontSize', FONT_SIZE_LABEL);
title('(c) Load rate distribution (KDE)', 'FontSize', FONT_SIZE_SUBPLOT);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

% Subplot 4: Performance radar chart
pax = polaraxes('Position', [0.55, 0.1, 0.4, 0.35]);
hold(pax, 'on');

% Calculate performance metrics
metrics_baseline = [
    mean(baseline_loads*100);
    100 - std(baseline_loads*100);
    100 - baseline_violate_num/baseline_NV*100;
    (100 - mean(baseline_route_lengths_local)/max(baseline_route_lengths_local)*100);
    100 - baseline_NV/cusnum*100;
];

if num_user_stations > 0
    metrics_hybrid = [
        mean(hybrid_loads*100);
        100 - std(hybrid_loads*100);
        100 - best_violate_num/bestNV*100;
        (100 - mean(hybrid_route_lengths_local)/max(hybrid_route_lengths_local)*100);
        100 - bestNV/cusnum*100;
    ];
else
    metrics_hybrid = metrics_baseline;
end

% Normalize
metrics_max = max([metrics_baseline, metrics_hybrid], [], 2);
metrics_baseline_norm = metrics_baseline ./ metrics_max * 100;
metrics_hybrid_norm = metrics_hybrid ./ metrics_max * 100;

theta = linspace(0, 2*pi, 6);
labels = {'Load Rate', 'Stability', 'Feasibility', 'Compactness', 'Vehicle Util.'};

polarplot(pax, theta, [metrics_baseline_norm; metrics_baseline_norm(1)], ...
          '-o', 'Color', COLOR_SCHEME.blue, 'LineWidth', 2, ...
          'MarkerSize', 8, 'MarkerFaceColor', COLOR_SCHEME.blue, ...
          'DisplayName', 'Baseline');

if num_user_stations > 0
    polarplot(pax, theta, [metrics_hybrid_norm; metrics_hybrid_norm(1)], ...
              '-^', 'Color', COLOR_SCHEME.green, 'LineWidth', 2, ...
              'MarkerSize', 8, 'MarkerFaceColor', COLOR_SCHEME.green, ...
              'DisplayName', 'Hybrid');
end

thetaticks(pax, rad2deg(theta(1:5)));
thetaticklabels(pax, labels);
title(pax, '(d) Performance radar chart', 'FontSize', FONT_SIZE_SUBPLOT);
legend(pax, 'Location', 'best', 'FontSize', FONT_SIZE_LEGEND);

fig10_path = fullfile(figures_dir, 'Fig10_3D_Efficiency_Analysis');
savefig(gcf, [fig10_path '.fig']);
fprintf('  Saved: Fig10_3D_Efficiency_Analysis.fig\n');

%% =========================================================================
% Figure 12: Algorithm Convergence Dynamics
% =========================================================================
fprintf('\nGenerating Figure 12: Convergence Dynamics...\n');
figure('Position', [250, 100, 1200, 500], 'Color', 'w', 'PaperPositionMode', 'auto');

% Subplot 1: Log-scale convergence curve
subplot(1, 2, 1);
hold on; box on;

semilogy(baseline_history.generation, baseline_history.best_obj, '-', ...
         'Color', COLOR_SCHEME.blue, 'LineWidth', 2, 'DisplayName', 'Baseline');

if num_user_stations > 0
    max_gen = min(length(baseline_history.generation), length(hybrid_history.generation));
    semilogy(hybrid_history.generation(1:max_gen), hybrid_history.best_obj(1:max_gen), '--', ...
             'Color', COLOR_SCHEME.green, 'LineWidth', 2, 'DisplayName', 'Hybrid');
end

xlabel('Generation', 'FontSize', FONT_SIZE_LABEL);
ylabel('Best objective (log scale)', 'FontSize', FONT_SIZE_LABEL);
title('(a) Convergence curve (log scale)', 'FontSize', FONT_SIZE_SUBPLOT);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

% Subplot 2: Population diversity evolution
subplot(1, 2, 2);
hold on; box on;

% Population diversity (gap between average and best)
diversity_baseline = (baseline_history.avg_obj - baseline_history.best_obj) ./ ...
                     baseline_history.best_obj * 100;

area(baseline_history.generation, diversity_baseline, 'FaceColor', COLOR_SCHEME.blue, ...
     'FaceAlpha', 0.3, 'EdgeColor', COLOR_SCHEME.blue, 'LineWidth', 2, ...
     'DisplayName', 'Baseline diversity');

if num_user_stations > 0
    diversity_hybrid = (hybrid_history.avg_obj(1:max_gen) - hybrid_history.best_obj(1:max_gen)) ./ ...
                       hybrid_history.best_obj(1:max_gen) * 100;
    area(hybrid_history.generation(1:max_gen), diversity_hybrid, 'FaceColor', COLOR_SCHEME.green, ...
         'FaceAlpha', 0.3, 'EdgeColor', COLOR_SCHEME.green, 'LineWidth', 2, ...
         'DisplayName', 'Hybrid diversity');
end

xlabel('Generation', 'FontSize', FONT_SIZE_LABEL);
ylabel('Population diversity (%)', 'FontSize', FONT_SIZE_LABEL);
title('(c) Population diversity evolution', 'FontSize', FONT_SIZE_SUBPLOT);
legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

%% =========================================================================
% Figure 13: Multi-dimensional Station Feature Analysis
% =========================================================================
fprintf('\nGenerating Figure 13: Station Features...\n');
figure('Position', [300, 100, 1200, 900], 'Color', 'w', 'PaperPositionMode', 'auto');

% Subplot 1: Bubble chart - distance vs demand
subplot(2, 2, 1);
hold on; box on;

if num_user_stations > 0
    % Operator stations
    operator_dist = dist_from_depot(operator_customers);
    operator_demand = demands(operator_customers);
    scatter(operator_dist, operator_demand, 80, COLOR_SCHEME.blue, 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.2, 'MarkerFaceAlpha', 0.6, ...
            'DisplayName', 'Operator stations');
    
    % User stations - bubble size = suitability score
    user_dist = dist_from_depot(selected_stations);
    user_demand = demands(selected_stations);
    user_scores = user_suitability_score(selected_stations);
    
    bubble_sizes = 50 + user_scores(:) .* 200;
    scatter(user_dist, user_demand, bubble_sizes, COLOR_SCHEME.green, 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.2, 'MarkerFaceAlpha', 0.6, ...
            'DisplayName', 'User stations');
    
    % Trend lines
    p_op = polyfit(operator_dist, operator_demand, 1);
    p_user = polyfit(user_dist, user_demand, 1);
    x_trend = linspace(min(dist_from_depot), max(dist_from_depot), 100);
    plot(x_trend, polyval(p_op, x_trend), '--', 'Color', COLOR_SCHEME.blue, ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(x_trend, polyval(p_user, x_trend), '--', 'Color', COLOR_SCHEME.green, ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    scatter(NaN, NaN, 150, 'k', 'DisplayName', 'Bubble size = Suitability');
    
    legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
else
    scatter(dist_from_depot, demands, 80, COLOR_SCHEME.blue, 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.2, 'MarkerFaceAlpha', 0.6);
end

xlabel('Distance from depot (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Demand (bikes)', 'FontSize', FONT_SIZE_LABEL);
title('(a) Distance-demand bubble chart', 'FontSize', FONT_SIZE_SUBPLOT);
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

% Subplot 2: Feature comparison (bar chart)
subplot(2, 2, 2);
hold on; box on;

if num_user_stations > 0
    feature_names = {'Distance\n(km)', 'Demand\n(bikes)', 'Abs Gap\n(bikes)', 'Incentive\n(CNY)'};
    
    % Mean values
    op_means = [
        mean(dist_from_depot(operator_customers)), ...
        mean(abs(demands(operator_customers))), ...
        mean(abs(supply_demand_gap(operator_customers))), ...
        0
    ];
    
    user_means = [
        mean(dist_from_depot(selected_stations)), ...
        mean(abs(demands(selected_stations))), ...
        mean(abs(supply_demand_gap(selected_stations))), ...
        mean(user_incentive_per_station(selected_stations))
    ];
    
    % Normalize
    max_vals = max([op_means; user_means], [], 1);
    op_norm = op_means ./ max_vals * 100;
    user_norm = user_means ./ max_vals * 100;
    
    % Grouped bar chart
    x = 1:4;
    width = 0.35;
    
    b1 = bar(x - width/2, op_norm, width, 'FaceColor', COLOR_SCHEME.blue, ...
             'EdgeColor', 'k', 'LineWidth', 1.2);
    b2 = bar(x + width/2, user_norm, width, 'FaceColor', COLOR_SCHEME.green, ...
             'EdgeColor', 'k', 'LineWidth', 1.2);
    
    % Value labels
    for i = 1:4
        text(i - width/2, op_norm(i) + 3, sprintf('%.1f', op_means(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        text(i + width/2, user_norm(i) + 3, sprintf('%.1f', user_means(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', ...
             'Color', [0 0.5 0]);
    end
    
    set(gca, 'XTick', x, 'XTickLabel', feature_names);
    ylim([0 110]);
    
    legend([b1, b2], {'Operator stations', 'User stations'}, ...
           'Location', 'northwest', 'FontSize', FONT_SIZE_LEGEND);
else
    text(0.5, 0.5, 'No user stations', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', FONT_SIZE_LABEL);
end

ylabel('Normalized value (%)', 'FontSize', FONT_SIZE_LABEL);
title('(b) Feature comparison (operator vs user stations)', 'FontSize', FONT_SIZE_SUBPLOT);
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

% Subplot 3: Decision boundary visualization
subplot(2, 2, 3);
hold on; box on;

if num_user_stations > 0
    % Grid
    dist_range = linspace(min(dist_from_depot), max(dist_from_depot), 50);
    demand_range = linspace(min(demands), max(demands), 50);
    [D_grid, Dem_grid] = meshgrid(dist_range, demand_range);
    
    % Decision function (user suitability)
    Z_decision = zeros(size(D_grid));
    for i = 1:size(D_grid, 1)
        for j = 1:size(D_grid, 2)
            remoteness = D_grid(i,j) / max(dist_from_depot);
            low_demand = 1 - Dem_grid(i,j) / max(demands);
            Z_decision(i,j) = 0.7 * remoteness + 0.3 * low_demand;
        end
    end
    
    % Decision boundary
    contourf(D_grid, Dem_grid, Z_decision, 20, 'LineStyle', 'none');
    colormap(flipud(gray));
    cb = colorbar;
    ylabel(cb, 'User suitability', 'FontSize', FONT_SIZE_LEGEND);
    
    % Overlay stations
    scatter(dist_from_depot(operator_customers), demands(operator_customers), ...
            60, COLOR_SCHEME.blue, 'o', 'filled', 'MarkerEdgeColor', 'w', ...
            'LineWidth', 1.5, 'DisplayName', 'Operator');
    scatter(dist_from_depot(selected_stations), demands(selected_stations), ...
            60, COLOR_SCHEME.green, '^', 'filled', 'MarkerEdgeColor', 'w', ...
            'LineWidth', 1.5, 'DisplayName', 'User');
    
    % Decision boundary line (threshold=0.5)
    contour(D_grid, Dem_grid, Z_decision, [0.5 0.5], 'r--', 'LineWidth', 2.5, ...
            'DisplayName', 'Decision boundary');
    
    legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
end

xlabel('Distance from depot (km)', 'FontSize', FONT_SIZE_LABEL);
ylabel('Demand (bikes)', 'FontSize', FONT_SIZE_LABEL);
title('(c) Dispatch decision boundary', 'FontSize', FONT_SIZE_SUBPLOT);

% Subplot 4: Feature correlation matrix heatmap
subplot(2, 2, 4);

if num_user_stations > 0
    % Feature correlation matrix
    feature_matrix = [dist_from_depot, demands, supply, supply_demand_gap];
    
    % Subset to avoid large computation
    subset_size = min(50, cusnum);
    subset_idx = round(linspace(1, cusnum, subset_size));
    feature_subset = feature_matrix(subset_idx, :);
    
    % Correlation matrix
    corr_matrix = corrcoef(feature_subset);
    
    imagesc(corr_matrix);
    colormap(jet);
    cb = colorbar;
    ylabel(cb, 'Correlation', 'FontSize', FONT_SIZE_LEGEND);
    caxis([-1, 1]);
    
    % Value labels
    for i = 1:4
        for j = 1:4
            text(j, i, sprintf('%.2f', corr_matrix(i,j)), ...
                 'HorizontalAlignment', 'center', 'FontSize', FONT_SIZE_LEGEND, ...
                 'FontWeight', 'bold', 'Color', 'w');
        end
    end
    
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Distance', 'Demand', 'Supply', 'Gap'});
    set(gca, 'YTick', 1:4, 'YTickLabel', {'Distance', 'Demand', 'Supply', 'Gap'});
end

title('(d) Feature correlation matrix', 'FontSize', FONT_SIZE_SUBPLOT);

sgtitle('Fig. 13 Multi-dimensional station characteristics and dispatch decision', ...
        'FontSize', FONT_SIZE_TITLE, 'FontWeight', 'normal');

fig13_path = fullfile(figures_dir, 'Fig13_Station_Features_Analysis');
savefig(gcf, [fig13_path '.fig']);
fprintf('  Saved: Fig13_Station_Features_Analysis.fig\n');

%% =========================================================================
% Figure 14: Route Length Distribution Comparison
% =========================================================================
fprintf('\nGenerating Figure 14: Route Length Distribution...\n');
figure('Position', [350, 100, 800, 500], 'Color', 'w', 'PaperPositionMode', 'auto');
hold on; box on;

% Calculate route lengths
baseline_route_lengths = zeros(baseline_NV, 1);
for i = 1:baseline_NV
    baseline_route_lengths(i) = part_length(baseline_VC{i}, dist);
end

if num_user_stations > 0
    hybrid_route_lengths = zeros(bestNV, 1);
    for i = 1:bestNV
        global_route = operator_customers(bestVC{i});
        hybrid_route_lengths(i) = part_length(global_route, dist);
    end
    
    % Comparison plot
    max_routes = max(baseline_NV, bestNV);
    x_baseline = 1:baseline_NV;
    x_hybrid = 1:bestNV;
    
    bar(x_baseline, baseline_route_lengths, 0.4, 'FaceColor', COLOR_SCHEME.blue, ...
        'EdgeColor', 'w', 'DisplayName', 'Baseline');
    bar(x_hybrid + 0.4, hybrid_route_lengths, 0.4, 'FaceColor', COLOR_SCHEME.green, ...
        'EdgeColor', 'w', 'DisplayName', 'Hybrid');
    
    legend('Location', 'best', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
else
    hybrid_route_lengths = baseline_route_lengths;
    bar(1:baseline_NV, baseline_route_lengths, 'FaceColor', COLOR_SCHEME.blue, 'EdgeColor', 'w');
end

xlabel('Route ID', 'FontSize', FONT_SIZE_LABEL);
ylabel('Route length (km)', 'FontSize', FONT_SIZE_LABEL);
title('Fig. 14 Route length distribution', 'FontSize', FONT_SIZE_TITLE, 'FontWeight', 'normal');
grid on;
set(gca, 'GridColor', COLOR_SCHEME.grid, 'GridAlpha', 0.3);

fig14_path = fullfile(figures_dir, 'Fig14_Route_Length_Distribution');
savefig(gcf, [fig14_path '.fig']);
fprintf('  Saved: Fig14_Route_Length_Distribution.fig\n');

%% =========================================================================
% Table Generation
% =========================================================================
fprintf('\nGenerating Tables...\n');

% Table 1: Baseline route details
fprintf('  Generating Table 1: Baseline Routes...\n');
table1_data = generate_route_table(baseline_VC, demands, dist, cap, a, b, L, s, ...
                                    fixed_cost_per_vehicle, variable_cost_per_km, ...
                                    service_cost_per_station, driver_hourly_cost, avg_speed);
writetable(table1_data, fullfile(tables_dir, 'Table1_Baseline_Routes.xlsx'));

% Table 2: Hybrid route details
fprintf('  Generating Table 2: Hybrid Routes...\n');
if num_user_stations > 0
    table2_data = generate_route_table(bestVC, demands(operator_customers), sub_dist, cap, ...
                                        a(operator_customers), b(operator_customers), L, ...
                                        s(operator_customers), fixed_cost_per_vehicle, ...
                                        variable_cost_per_km, service_cost_per_station, ...
                                        driver_hourly_cost, avg_speed, operator_customers);
else
    table2_data = table1_data;
end
writetable(table2_data, fullfile(tables_dir, 'Table2_Hybrid_Routes.xlsx'));

% Table 3: User station information
fprintf('  Generating Table 3: User Stations...\n');
if num_user_stations > 0
    table3_data = table();
    table3_data.StationID = selected_stations;
    table3_data.X_Coordinate = customer(selected_stations, 1);
    table3_data.Y_Coordinate = customer(selected_stations, 2);
    table3_data.Demand = demands(selected_stations);
    table3_data.Supply = supply(selected_stations);
    table3_data.Inventory = supply(selected_stations);
    table3_data.DemandLevel = supply_demand_gap(selected_stations);
    table3_data.IncentiveCost = user_incentive_per_station(selected_stations);
    table3_data.SuitabilityScore = user_suitability_score(selected_stations);
else
    table3_data = table({'No user stations'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, ...
                        'VariableNames', {'StationID', 'X_Coordinate', 'Y_Coordinate', ...
                        'Demand', 'Supply', 'Inventory', 'DemandLevel', ...
                        'IncentiveCost', 'SuitabilityScore'});
end
writetable(table3_data, fullfile(tables_dir, 'Table3_User_Stations.xlsx'));

% Table 5: Overall comparison
fprintf('  Generating Table 5: Overall Comparison...\n');
table5_data = table();
table5_data.Method = {'Baseline (Operator-only)'; 'Hybrid (Operator+User)'};
table5_data.NumTrucks = [baseline_NV; bestNV];
table5_data.RouteDistance_km = [baseline_TD; bestTD];
table5_data.NumStationsServiced = [cusnum; cusnum - num_user_stations];
table5_data.TotalCost_CNY = [baseline_total_cost; final_total_cost];
avg_baseline_load = mean(baseline_loads(~isnan(baseline_loads)));
avg_hybrid_load = mean(hybrid_loads(~isnan(hybrid_loads)));
table5_data.AvgLoadRate_percent = [avg_baseline_load * 100; avg_hybrid_load * 100];
table5_data.Optimization_percent = [0; (baseline_total_cost - final_total_cost) / baseline_total_cost * 100];
writetable(table5_data, fullfile(tables_dir, 'Table5_Overall_Comparison.xlsx'));

%% =========================================================================
% Save Workspace
% =========================================================================
save('my_workspace.mat');
fprintf('\n=================================================================\n');
fprintf('All processing complete!\n');
fprintf('=================================================================\n');
toc