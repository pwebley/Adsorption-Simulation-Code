
function [nodes, valves] = define_valves(tspan)
% DEFINE_VALVES   Defines the flow network and valve configurations

% === Define nodes with pressure (Pa), temperature (K), and composition ===
nodes.FeedTank.P     = 1e5;         % Feed pressure (Pa)
nodes.FeedTank.T     = 318;         % Feed temperature [K]
nodes.FeedTank.y     = [0.003, 0.21, 0.787];  % Mole fractions of (H2O, O2, N2)
nodes.ProductTank.P  = 9e4;         %ProductTank Pressure (Pa)
nodes.ProductTank.T  = 318;         %Product Tank T (K)
nodes.ProductTank.y  = [0.00,1.00,0.00]; %Product tank composition
nodes.Waste.P        = 1e4;

% Define valves and assign Cv value when open
valve_ids = {'V1','V2','V3','V4','V5','V6','V7','V8','V9'};
descriptions = {
    'Feed to Bed1', 'Feed to Bed2', 'Waste from Bed1', 'Waste from Bed2', ...
    'Bed1 to Product', 'Bed2 to Product', 'Bed 1 to/from Bed2','Purge to Bed1',...
    'Purge to Bed2'};
from_nodes = {
    'FeedTank', 'FeedTank', 'Bed1_inlet', 'Bed2_inlet', ...
     'Bed1_outlet', 'Bed2_outlet', 'Bed1_outlet','ProductTank',...
     'ProductTank'
     };
to_nodes = {
    'Bed1_inlet', 'Bed2_inlet', 'Waste', 'Waste', ...
    'ProductTank', 'ProductTank', 'Bed2_outlet','Bed1_outlet',...
    'Bed2_outlet'
};

Cv_vals = [NaN, NaN, 0.015, 0.015, 0.0008, 0.0008, 0.0012,0.0001,0.0001];  % NaN = flow controlled

for i = 1:length(valve_ids)
    valves(i).id = valve_ids{i};
    valves(i).from = from_nodes{i};
    valves(i).to = to_nodes{i};
    valves(i).allow_reverse = false;
    valves(i).description = descriptions{i};

    if isnan(Cv_vals(i))
        valves(i).mode = 'flow_controlled';
        valves(i).Cv_max = 0;
        valves(i).fixed_flow = 0.1;  % adjust if needed
    else
        valves(i).mode = 'Cv';
        valves(i).Cv_max = Cv_vals(i);
        valves(i).fixed_flow = 0;
    end
end
% Designate the PE valve as reverse flow
    valves(7).allow_reverse = true;

end
