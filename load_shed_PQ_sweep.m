function [] = load_shed_PQ_sweep(casedata, PQ_sweep_data)

if nargin < 2
    PQ_sweep_data = sprintf('test-convergence-%s.csv', casedata);
end

define_constants;
[o_BUS_I, o_S, o_theta, o_P, o_Q, o_VM, o_Pi, o_Qi, o_PD_tot, o_QD_tot, o_PG_orig, o_QG_orig] = deal(1,2,3,4,5,6,7,8,9,10,11,12);

if nargin < 2
    PQ_sweep_data = sprintf('test-convergence-%s.csv', casedata);
end

filename = sprintf('load-shed-PQ-sweep-%s.csv', casedata);
fileID = fopen(filename, 'w');

PQ_sweep = csvread(PQ_sweep_data, 1);
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%si\n','Bus_I', 'Bus_type', 'min_bus', 'min_S', 'min_ang', 'num_soln', 'num_attempts')
mpc = loadcase(casedata);

num_attacks = size(PQ_sweep, 1);

for a = 1:num_attacks
    attack_mpc = loadcase(casedata);
    P = PQ_sweep(a, o_P);
    Q = PQ_sweep(a, o_Q);
    bus_idx = find(attack_mpc.bus(:, BUS_I) == PQ_sweep(a, o_BUS_I));
    attack_mpc.bus(bus_idx, [PD, QD]) = [P, Q];
    num_soln = 0;
    num_buses = size(mpc.bus, 1);
    min_S = inf;
    min_bus = 0;
    min_ang = 0;
    num_attempts = 0;
    for b = 1:num_buses
        shed_mpc = attack_mpc;
        bus_P = shed_mpc.bus(b, PD);
        bus_Q = shed_mpc.bus(b, QD);
        bus_S = sqrt(bus_P^2 + bus_Q^2);
        bus_ang = atan(bus_Q / bus_P);
        if bus_S ~= 0
            num_attempts = num_attempts + 1;
            fprintf('Attack %i / %i: Shed %i / %i\n', a, num_attacks, b, size(shed_mpc.bus, 1));
            shed_mpc.bus(b, [PD, QD]) = [0, 0];
            results = runpf(shed_mpc, mpoption('verbose', 0, 'out.all', 0));
            if results.success == 0
                if bus_S < min_S
                    min_S = bus_S;
                    min_bus = b;
                    min_ang = bus_ang;
                num_soln = num_soln + 1;
                end
            end
        end
    end
    fprintf(fileID, '%i,%i,%i,%i,%i,%i,%i\n', mpc.bus(bus_idx, BUS_I), 1, min_bus, min_S, min_ang, num_soln, num_attempts); 
end
fclose(fileID); 
                

