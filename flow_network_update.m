function [nodes, Q_valves, inlet_bc, outlet_bc] = flow_network_update(t, y, parm_all, valves, nodes, sim)
% FLOW_NETWORK_UPDATE   Updates node pressures and computes valve flows

% === Extract interface states ===
[P1_in, P1_out, T1_in, T1_out, y1_in, y1_out] = get_bed_interface_props(y, 1, parm_all{1});
[P2_in, P2_out, T2_in, T2_out, y2_in, y2_out] = get_bed_interface_props(y, 2, parm_all{2});
R = parm_all{1}.R;  % Gas constant

nodes.Bed1_inlet  = struct('P', P1_in,  'T', T1_in,  'y', y1_in);
nodes.Bed1_outlet = struct('P', P1_out, 'T', T1_out, 'y', y1_out);
nodes.Bed2_inlet  = struct('P', P2_in,  'T', T2_in,  'y', y2_in);
nodes.Bed2_outlet = struct('P', P2_out, 'T', T2_out, 'y', y2_out);

% === Initialize outputs ===
Q_valves = repmat(struct('id', '', 'Q', 0), length(valves), 1);
inlet_bc = cell(1, 2);
outlet_bc = cell(1, 2);

% === Default: zero-flow BCs using current node properties ===
inlet_bc{1}  = struct('Q', 0, 'T', T1_in,  'y', y1_in,  'Ct', P1_in / (R * T1_in));
outlet_bc{1} = struct('Q', 0, 'T', T1_out, 'y', y1_out, 'Ct', P1_out / (R * T1_out));
inlet_bc{2}  = struct('Q', 0, 'T', T2_in,  'y', y2_in,  'Ct', P2_in / (R * T2_in));
outlet_bc{2} = struct('Q', 0, 'T', T2_out, 'y', y2_out, 'Ct', P2_out / (R * T2_out));

% === Compute smoothed Cv and flow through each valve ===
for v = 1:length(valves)
    valve = valves(v);

    % Extract Cv schedule and smoothing info
    Cv_vec = valve.Cv_schedule;
    interval_start = valve.interval_start;
    interval_end = valve.interval_end;
    smoothing_duration = valve.smoothing_duration;

    % Compute smoothed Cv at time t using logistic transition
    Cv_t = 0;
    for i = 1:length(Cv_vec)
        t0 = double(interval_start(i));
        t1 = double(interval_end(i));
        if t >= t0 && t < t1
            if t < t0 + smoothing_duration
                % Smooth transition
                k = 4 / (smoothing_duration + eps);
                center = t0 + smoothing_duration / 2;
                weight = 1 / (1 + exp(-k * (t - center)));
                Cv_t = weight * Cv_vec(i);
            else
                % Fully open
                Cv_t = Cv_vec(i);
            end
            break
        end
    end

    if Cv_t < 1e-12
        continue  % Valve is closed → skip flow and BC logic
    end

    % Compute pressure drop and flow
    Pin = nodes.(valve.from).P;
    Pout = nodes.(valve.to).P;
    dP = Pin - Pout;
    Q = Cv_t * sqrt(abs(dP)) * sign(dP);

    if ~valve.allow_reverse && Q < 0
        Q = 0;
    end

    Q_valves(v).id = valve.id;
    Q_valves(v).Q = Q;

    % === Handle inter-bed outlet ↔ outlet connections ===
    if (strcmp(valve.from, 'Bed1_outlet') && strcmp(valve.to, 'Bed2_outlet')) || ...
       (strcmp(valve.from, 'Bed2_outlet') && strcmp(valve.to, 'Bed1_outlet'))

        if Q >= 0
            outlet_bc{1}.Q =  abs(Q);
            outlet_bc{2}.Q = -abs(Q);
            outlet_bc{1}.T = T1_out;
            outlet_bc{1}.y = y1_out;
            outlet_bc{1}.Ct = P1_out / (R * T1_out);
            outlet_bc{2}.T = T1_out;
            outlet_bc{2}.y = y1_out;
            outlet_bc{2}.Ct = P1_out / (R * T1_out);
        else
            outlet_bc{1}.Q = -abs(Q);
            outlet_bc{2}.Q =  abs(Q);
            outlet_bc{1}.T = T2_out;
            outlet_bc{1}.y = y2_out;
            outlet_bc{1}.Ct = P2_out / (R * T2_out);
            outlet_bc{2}.T = T2_out;
            outlet_bc{2}.y = y2_out;
            outlet_bc{2}.Ct = P2_out / (R * T2_out);
        end
        continue
    end

    % === Generic boundary condition handling ===
    T_valve = nodes.(valve.from).T;
    y_valve = nodes.(valve.from).y;
    P_valve = nodes.(valve.from).P;

    if strcmp(valve.to, 'Bed1_inlet')
        inlet_bc{1}.Q = abs(Q);
        inlet_bc{1}.T = T_valve;
        inlet_bc{1}.y = y_valve;
        inlet_bc{1}.Ct = P_valve / (R * T_valve);
    elseif strcmp(valve.from, 'Bed1_inlet')
        inlet_bc{1}.Q = -abs(Q);
        inlet_bc{1}.T = T_valve;
        inlet_bc{1}.y = y_valve;
        inlet_bc{1}.Ct = P_valve / (R * T_valve);
    end

    if strcmp(valve.to, 'Bed2_inlet')
        inlet_bc{2}.Q = abs(Q);
        inlet_bc{2}.T = T_valve;
        inlet_bc{2}.y = y_valve;
        inlet_bc{2}.Ct = P_valve / (R * T_valve);
    elseif strcmp(valve.from, 'Bed2_inlet')
        inlet_bc{2}.Q = -abs(Q);
        inlet_bc{2}.T = T_valve;
        inlet_bc{2}.y = y_valve;
        inlet_bc{2}.Ct = P_valve / (R * T_valve);
    end

    if strcmp(valve.from, 'Bed1_outlet')
        outlet_bc{1}.Q = abs(Q);
        outlet_bc{1}.T = T_valve;
        outlet_bc{1}.y = y_valve;
        outlet_bc{1}.Ct = P_valve / (R * T_valve);
    elseif strcmp(valve.to, 'Bed1_outlet')
        outlet_bc{1}.Q = -abs(Q);
        outlet_bc{1}.T = T_valve;
        outlet_bc{1}.y = y_valve;
        outlet_bc{1}.Ct = P_valve / (R * T_valve);
    end

    if strcmp(valve.from, 'Bed2_outlet')
        outlet_bc{2}.Q = abs(Q);
        outlet_bc{2}.T = T_valve;
        outlet_bc{2}.y = y_valve;
        outlet_bc{2}.Ct = P_valve / (R * T_valve);
    elseif strcmp(valve.to, 'Bed2_outlet')
        outlet_bc{2}.Q = -abs(Q);
        outlet_bc{2}.T = T_valve;
        outlet_bc{2}.y = y_valve;
        outlet_bc{2}.Ct = P_valve / (R * T_valve);
    end
end
end
