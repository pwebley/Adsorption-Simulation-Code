
function valves = define_valves()
% DEFINE_VALVES   Defines the flow network and valve configurations

% === Define nodes with pressure (Pa), temperature (K), and composition ===
nodes.FeedTank.P     = 1e5;         % Feed pressure (Pa)
nodes.FeedTank.T     = 318;         % Feed temperature [K]
nodes.FeedTank.y     = [0.003, 0.21, 0.787];  % Mole fractions of (H2O, O2, N2)
nodes.ProductTank.P  = 9e4;         % Product Tank Pressure (Pa)
nodes.ProductTank.T  = 318;         % Product Tank Temperature (K)
nodes.ProductTank.y  = [0.00, 1.00, 0.00]; % Product tank composition
nodes.Waste.P        = 1e4;

% === Define valves and connectivity ===
valve_ids = {'V1','V2','V3','V4','V5','V6','V7','V8','V9'};
descriptions = {
    'Feed to Bed1', 'Feed to Bed2', 'Waste from Bed1', 'Waste from Bed2', ...
    'Bed1 to Product', 'Bed2 to Product', 'Bed 1 to/from Bed2','Purge to Bed1',...
    'Purge to Bed2'};
from_nodes = {
    'FeedTank', 'FeedTank', 'Bed1_inlet', 'Bed2_inlet', ...
     'Bed1_outlet', 'Bed2_outlet', 'Bed1_outlet','ProductTank',...
     'ProductTank'};
to_nodes = {
    'Bed1_inlet', 'Bed2_inlet', 'Waste', 'Waste', ...
    'ProductTank', 'ProductTank', 'Bed2_outlet','Bed1_outlet',...
    'Bed2_outlet'};

% === Initialize valve structs ===
for i = 1:length(valve_ids)
    valves(i).id = valve_ids{i};
    valves(i).from = from_nodes{i};
    valves(i).to = to_nodes{i};
    valves(i).allow_reverse = false;
    valves(i).description = descriptions{i};
    valves(i).mode = 'Cv_schedule';   % Using time-varying Cv schedule
end

% Mark valve V7 (index 7) as reversible
valves(7).allow_reverse = true;

% === Load Cv matrix and time intervals ===
data = load('valve_cv_matrix_data.mat');
Cv_matrix = data.Cv_matrix;              % [n_intervals x n_valves]
interval_start = data.interval_start(:); % column vector
interval_end = data.interval_end(:);     % column vector

% Assign Cv schedule and shared timing to each valve
for i = 1:length(valves)
    valves(i).Cv_schedule = Cv_matrix(:, i+1);
    valves(i).interval_start = interval_start;
    valves(i).interval_end = interval_end;
    valves(i).smoothing_duration = 1.0;  % seconds for logistic smoothing
end

% === Pass node definitions out via valve struct ===
valves(1).nodes = nodes;  % included only once to avoid redundancy
end
