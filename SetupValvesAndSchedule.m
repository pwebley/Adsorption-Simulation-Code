function [nodes, valves, step_times, valve_schedule] = SetupValvesAndSchedule()
%% === Node Network Definition ===
sim.nodes = struct();

sim.nodes(1).name = 'exit1';
sim.nodes(1).type = 'fixed_pressure';
sim.nodes(1).pressure = 101325;
sim.nodes(1).temperature = 298.15;

sim.nodes(2).name = 'exit2';
sim.nodes(2).type = 'fixed_pressure';
sim.nodes(2).pressure = 101325;
sim.nodes(2).temperature = 298.15;

%% === Valve and node Network Definition ===
sim.num_valves = 20;

% Define each valve with connectivity and mode
valves(1).name = 'V1';
valves(1).from = 'Tank1';
valves(1).to = 'Bed1_in';
valves(1).mode = 'fixed_flow';  % 'Cv' | 'fixed_flow' | 'vacuum_pump' | 'compressor'
valves(1).fixed_flow = 0.005;  % mol/s
valves(1).Cv = NaN;
valves(1).reverse = false;

valves(2).name = 'V2';
valves(2).from = 'Tank1';
valves(2).to = 'Bed2_in';
valves(2).mode = 'fixed_flow';  % 'Cv' | 'fixed_flow' | 'vacuum_pump' | 'compressor'
valves(2).fixed_flow = 0.005;  % mol/s
valves(2).Cv = NaN;
valves(2).reverse = false;

valves(3).name = 'V3';
valves(3).from = 'Tank1';
valves(3).to = 'Bed3_in';
valves(3).mode = 'fixed_flow';  % 'Cv' | 'fixed_flow' | 'vacuum_pump' | 'compressor'
valves(3).fixed_flow = 0.005;  % mol/s
valves(3).Cv = NaN;
valves(3).reverse = false;

valves(4).name = 'V4';
valves(4).from = 'Tank2';
valves(4).to = 'Bed1_in';
valves(4).mode = 'fixed_flow';  % 'Cv' | 'fixed_flow' | 'vacuum_pump' | 'compressor'
valves(4).fixed_flow = 0.001;  % mol/s
valves(4).Cv = NaN;
valves(4).reverse = false;

valves(5).name = 'V5';
valves(5).from = 'Tank2';
valves(5).to = 'Bed2_in';
valves(5).mode = 'fixed_flow';  % 'Cv' | 'fixed_flow' | 'vacuum_pump' | 'compressor'
valves(5).fixed_flow = 0.001;  % mol/s
valves(5).Cv = NaN;
valves(5).reverse = false;

valves(6).name = 'V6';
valves(6).from = 'Tank2';
valves(6).to = 'Bed3_in';
valves(6).mode = 'fixed_flow';  % 'Cv' | 'fixed_flow' | 'vacuum_pump' | 'compressor'
valves(6).fixed_flow = 0.001;  % mol/s
valves(6).Cv = NaN;
valves(6).reverse = false;

valves(7).name = 'V7';
valves(7).from = 'Bed1_out';
valves(7).to = 'Tank2';
valves(7).mode = 'vacuum_pump';
valves(7).Cv = NaN;  % Not used for vacuum pump
valves(7).reverse = false;

valves(8).name = 'V8';
valves(8).from = 'Bed2_out';
valves(8).to = 'Tank2';
valves(8).mode = 'vacuum_pump';
valves(8).Cv = NaN;  % Not used for vacuum pump
valves(8).reverse = false;

valves(9).name = 'V9';
valves(9).from = 'Bed3_out';
valves(9).to = 'Tank2';
valves(9).mode = 'vacuum_pump';
valves(9).Cv = NaN;  % Not used for vacuum pump
valves(9).reverse = false;

valves(10).name = 'V10';
valves(10).from = 'Bed1_out';
valves(10).to = 'Tank3';
valves(10).mode = 'Cv';
valves(10).Cv = 0.001;  
valves(10).reverse = false;

valves(11).name = 'V11';
valves(11).from = 'Bed2_out';
valves(11).to = 'Tank3';
valves(11).mode = 'Cv';
valves(11).Cv = 0.001;  
valves(11).reverse = false;

valves(12).name = 'V12';
valves(12).from = 'Bed3_out';
valves(12).to = 'Tank3';
valves(12).mode = 'Cv';
valves(12).Cv = 0.001;  
valves(12).reverse = false;

valves(13).name = 'V13';
valves(13).from = 'Tank3';
valves(13).to = 'Bed1_out';
valves(13).mode = 'Cv';
valves(13).Cv = 0.001;  
valves(13).reverse = false;

valves(14).name = 'V14';
valves(14).from = 'Tank3';
valves(14).to = 'Bed2_out';
valves(14).mode = 'Cv';
valves(14).Cv = 0.001;  
valves(14).reverse = false;

valves(15).name = 'V15';
valves(15).from = 'Tank3';
valves(15).to = 'Bed3_out';
valves(15).mode = 'Cv';
valves(15).Cv = 0.001;  
valves(15).reverse = false;

valves(16).name = 'V16';
valves(16).from = 'Bed1_out';
valves(16).to = 'Bed2_out';
valves(16).mode = 'Cv';
valves(16).Cv = 0.0001;  
valves(16).reverse = true;

valves(17).name = 'V17';
valves(17).from = 'Bed2_out';
valves(17).to = 'Bed3_out';
valves(17).mode = 'Cv';
valves(17).Cv = 0.0001;  
valves(17).reverse = true;

valves(18).name = 'V18';
valves(18).from = 'Bed1_out';
valves(18).to = 'Bed3_out';
valves(18).mode = 'Cv';
valves(18).Cv = 0.0001;  
valves(18).reverse = true;

valves(19).name = 'V19';
valves(19).from = 'Tank2';
valves(19).to = 'exit1';
valves(19).mode = 'Cv';
valves(19).Cv = 0.0001;  
valves(19).reverse = false;

valves(20).name = 'V20';
valves(20).from = 'Tank3';
valves(20).to = 'exit2';
valves(20).mode = 'Cv';
valves(20).Cv = 0.0001;  
valves(20).reverse = false;

sim.valves = valves;
node_names = unique([string({valves.from}), string({valves.to})]);
sim.node_names = node_names;

%% === Valve Schedule (time-based step logic) ===
% Define time steps and valve status matrix (1 = open, 0 = closed)
sim.step_times = [0, 25, 45, 50, 55, 80, 100, 105, 110, 135, 155, 160, 165];  % [s] boundaries between steps (n_steps+1)
n_steps = numel(sim.step_times) - 1;

% Matrix: [n_steps x n_valves] of logicals
sim.valve_open = [
    1 1 0 0 0 0 0 0 0 0 0 0 ; %V1
    0 0 0 0 0 0 0 0 1 1 0 0 ; %V2
    0 0 0 0 1 1 0 0 0 0 0 0 ; %V3
    0 0 0 0 1 1 1 0 0 0 0 0 ; %V4
    1 1 1 0 0 0 0 0 0 0 0 0 ; %V5
    0 0 0 0 0 0 0 0 1 1 1 0 ; %V6
    0 0 0 0 0 0 0 0 1 1 1 0 ; %V7
    0 0 0 0 1 1 1 0 0 0 0 0 ; %V8
    1 1 1 0 0 0 0 0 0 0 0 0 ; %V9
    0 1 1 0 0 1 1 0 0 0 0 0 ; %V10
    0 1 1 0 0 0 0 0 0 1 1 0 ; %V11
    0 0 0 0 0 1 1 0 0 1 1 0 ; %V12
    0 0 0 0 0 0 0 0 0 0 0 0 ; %V13
    0 0 0 0 0 0 0 0 0 0 0 0 ; %V14
    0 0 0 0 0 0 0 0 0 0 0 0 ; %V15
    0 0 0 0 0 0 0 0 0 0 0 1 ; %V16
    0 0 0 0 0 0 0 1 0 0 0 0 ; %V17
    0 0 0 1 0 0 0 0 0 0 0 0 ; %V18
    1 1 1 1 1 1 1 1 1 1 1 1 ; %V19
    1 1 1 1 1 1 1 1 1 1 1 1 ; %V20
]';

% Preallocate a struct array
for i = 1:n_steps
    for j = 1:sim.num_valves
        sim.valve_schedule(i,j).active = sim.valve_open(i,j);
        sim.valve_schedule(i,j).value = NaN;  % default: use sim.valves(j).Cv or .fixed_flow
    end
end
% Set values of valve CV's, fixed_flow, vacuum flow, etc for each valve in
% each step. Only enter values for steps where the valve is open.

% V1 fixed_flow in step 1, 2
sim.valve_schedule(1, 1).value = 0.007;  % mol/s for fixed_flow
sim.valve_schedule(2, 1).value = 0.01;   % mol/s for fixed_flow

% V2 fixed_flow in step 9, 10
sim.valve_schedule(2, 2).value = 0.007;  % mol/s for fixed_flow
sim.valve_schedule(10, 2).value = 0.01;   % mol/s for fixed_flow

% V3 fixed_flow in step 5, 6
sim.valve_schedule(5, 3).value = 0.007;    % mol/s for fixed_flow
sim.valve_schedule(6, 3).value = 0.01;     % mol/s for fixed_flow

% V4 fixed_flow in step 5, 6, 7
sim.valve_schedule(5, 4).value = 0.002;    % mol/s for fixed_flow
sim.valve_schedule(6, 4).value = 0.002;    % mol/s for fixed_flow
sim.valve_schedule(7, 4).value = 0.002;    % mol/s for fixed_flow

% V5 fixed_flow in step 1, 2, 3
sim.valve_schedule(1, 5).value = 0.002;    % mol/s for fixed_flow
sim.valve_schedule(2, 5).value = 0.002;    % mol/s for fixed_flow
sim.valve_schedule(3, 5).value = 0.002;    % mol/s for fixed_flow

% V6 fixed_flow in step 9, 10, 11
sim.valve_schedule(9, 6).value = 0.002;    % mol/s for fixed_flow
sim.valve_schedule(10, 6).value = 0.002;   % mol/s for fixed_flow
sim.valve_schedule(11, 6).value = 0.002;   % mol/s for fixed_flow

% V7 vacuum flow in step 9, 10, 11
sim.valve_schedule(9, 7).value = 0.02;     % m3/s for vacuum flow
sim.valve_schedule(10, 7).value = 0.02;    % m3/s for vacuum flow
sim.valve_schedule(11, 7).value = 0.02;    % m3/s for vacuum flow

% V8 vacuum flow in step 5, 6, 7
sim.valve_schedule(5, 8).value = 0.02;     % m3/s for vacuum flow
sim.valve_schedule(6, 8).value = 0.02;     % m3/s for vacuum flow
sim.valve_schedule(7, 8).value = 0.02;     % m3/s for vacuum flow

% V9 vacuum flow in step 1, 2, 3
sim.valve_schedule(1, 9).value = 0.02;     % m3/s for vacuum flow
sim.valve_schedule(2, 9).value = 0.02;     % m3/s for vacuum flow
sim.valve_schedule(3, 9).value = 0.02;     % m3/s for vacuum flow

% V10 Cv flow in step 2, 3, 6, 7
sim.valve_schedule(2, 10).value = 0.02;    % Cv
sim.valve_schedule(3, 10).value = 0.002;   % Cv
sim.valve_schedule(6, 10).value = 0.005;   % Cv
sim.valve_schedule(7, 10).value = 0.005;   % Cv

% V11 Cv flow in step 2, 3, 10, 11
sim.valve_schedule(2, 11).value = 0.005;    % Cv
sim.valve_schedule(3, 11).value = 0.005;   % Cv
sim.valve_schedule(10, 11).value = 0.02;  % Cv
sim.valve_schedule(11, 11).value = 0.002;  % Cv

% V12 Cv flow in step 6, 7, 10, 11
sim.valve_schedule(6, 12).value = 0.02;    % Cv
sim.valve_schedule(7, 12).value = 0.002;   % Cv
sim.valve_schedule(10, 12).value = 0.005;  % Cv
sim.valve_schedule(11, 12).value = 0.005;  % Cv

% V13 to V15 not used in this cycle

% V16 CV flow bed1<->Bed2 Pressure equalisation, step 4
sim.valve_schedule(4, 16).value = 0.002;    % Cv

% V17 CV flow bed2<->Bed3 Pressure equalisation, step 8
sim.valve_schedule(8, 17).value = 0.002;    % Cv

% V18 CV flow bed1<->Bed3 Pressure equalisation, step 12
sim.valve_schedule(12, 18).value = 0.002;    % Cv

end
