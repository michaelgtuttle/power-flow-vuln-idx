function [out_file] = test_convergence(case_data, max_power, angular_res, radial_step)
reset(gca);
reset(gcf);
if nargin < 4
    radial_step = 10;
    if nargin < 3
        angular_res = 12;
        if nargin < 2
            max_power = inf;
        end
    end
end

opt = mpoption('verbose', 0, 'out.all', 0);

define_constants;
[o_BUS_I, o_S, o_theta, o_P, o_Q, o_VM, o_Pi, o_Qi, o_PD_tot, o_QD_tot, o_PG_orig, o_QG_orig] = deal(1,2,3,4,5,6,7,8,9,10,11,12);



mpc = loadcase(case_data);
orig_res = runpf(mpc, opt);

PG_orig = 0;
QG_orig = 0;
for i = 1:size(orig_res.gen, 1)
    PG_orig = PG_orig + orig_res.gen(i, PG);
    QG_orig = QG_orig + orig_res.gen(i, QG);
end

PD_orig = 0;
QD_orig = 0;
for i = 1:size(mpc.bus, 1)
    PD_orig = mpc.bus(i, PD) + PD_orig;
    QD_orig = QD_orig + mpc.bus(i, QD);
end


pq = find(mpc.bus(:,BUS_TYPE) == 1);

file_name = sprintf('test-convergence-%s.csv', case_data);

fileID = fopen(file_name, 'w');

S_mag = 0;
theta = 0;
theta_inc = 2 * pi / angular_res;
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'BUS_I', 'S', 'theta', 'P', 'Q', 'VM', 'P_i', 'Q_i', 'PD_tot', 'QD_tot','PG_orig','QG_orig');

num_pq = size(pq, 1);

Success = 0;

for i = 1:num_pq 
    S_mag = 0; 
    theta = 0;
    failed = 0;
    bus_idx = pq(i); 
    P_i = mpc.bus(bus_idx, PD);
    Q_i = mpc.bus(bus_idx, QD);
    S_i = sqrt(P_i ^ 2 + Q_i^2);
    results = orig_res; 
    while S_mag < max_power
        fprintf('BUS %i / %i: S = %ie^j%i\n', i, num_pq, S_mag, theta);
        new_mpc = mpc;
        S = S_mag * exp(1j * theta);
        P = real(S) + P_i;
        Q = imag(S) + Q_i;
        %S_f = sqrt(P^2 + Q^2);
        PD_fin = PD_orig - P_i + P;
        QD_fin = QD_orig - Q_i + Q;
        new_mpc.bus(bus_idx, PD) = P; 
        new_mpc.bus(bus_idx, QD) = Q;
        results = runpf(new_mpc, opt);
        if results.success == 0
            failed = 1;
            Success = 1;
            %BUS, S_mag, theta, P_f, Q_f, VM, P_i, Q_i, PDtot_f, QDtot_f, PGtot_orig, QGtot_orig
            fprintf(fileID, '%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n', mpc.bus(bus_idx, BUS_I), S_mag, theta, P, Q, results.bus(bus_idx, VM), P_i, Q_i, PD_fin, QD_fin, PG_orig, QG_orig); 
        end
        theta = mod(theta + theta_inc, 2 * pi);
        if abs(theta) < 0.0001 || abs(theta - 2 * pi) < 0.0001
             theta = 0;
             if failed == 1
                 break;
             end
             S_mag = S_mag + radial_step;
        end
    end
end
fclose(fileID);

if Success == 0
    error('No buses failed to converge')
end

mpc = loadcase(case_data);
out = csvread(file_name, 1);
powers = zeros(size(pq, 1), size(out, 2));
for i = 1:size(pq, 1)
    bus_idx = pq(i);
    bus_name = mpc.bus(bus_idx, BUS_I);
    bus_rows = find(out(:, 1) == bus_name);
    powers(i, :) = out(bus_rows(1), :);
end

fig = figure();
subplot(2,1,1);
line = linspace(0, 10, 11);
%histogram(powers(:,o_S) / mpc.baseMVA,'NumBins', 10, 'FaceColor', [1, .5, .5])
h = histogram(powers(:,o_S) / mpc.baseMVA,'BinEdges', line, 'FaceColor', [1, .5, .5])
histogram('BinEdges', h.BinEdges, 'BinCounts', h.BinCounts / size(powers, 1) * 100, 'FaceColor', [1, .5, .5])
xlabel('|S| (pu)');
ylabel('Percent Failed (%)');
set(gca, 'fontsize', 12)
subplot(2,1,2);
h = histogram(out(:,o_theta) / pi, 'NumBins', angular_res, 'FaceColor', [1 .5 .5])
histogram('BinEdges', h.BinEdges, 'BinCounts', h.BinCounts / size(out, 1) * 100, 'FaceColor', [1, .5, .5])
xlabel('\theta (rad * \pi)');
ylabel('Percent Failed (%)'); 
set(gca, 'fontsize', 12)
fig2 = figure();
VM_edges = linspace(0, 2, 20);
h = histogram(out(:, o_VM), 'BinEdges', VM_edges, 'FaceColor', [1, .5, .5]);
histogram('BinEdges', h.BinEdges, 'BinCounts', h.BinCounts / size(out, 1) * 100, 'FaceColor', [1, .5, .5])
xlabel('Vm');
ylabel('Percent Failed (%)');
figurename = sprintf('test-convergence-%s-hist', case_data);
figurename2 = sprintf('test-convergence-%s-VM-hist', case_data);
print(fig, figurename, '-dpdf');
print(fig2, figurename2, '-dpdf');
out_file = file_name;
end






