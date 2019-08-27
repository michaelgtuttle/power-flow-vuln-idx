function [] = slacktest_load_shed(casedata, slack_data)
if nargin < 2
    slack_data = sprintf('slacktest-full-%s-pf.csv', casedata);
end
opt = mpoption('verbose', 0, 'out.all', 0);
problem = @(x) runpf(x, opt);
filename = sprintf('load-shed-slacktest-%s.csv', casedata);
fileID = fopen(filename, 'w');
define_constants;
[s_BUS_I, s_BUS_TYPE, s_PD, s_QD, s_PG, s_QG, s_VM_BF, s_VM, s_VA_BPF, s_VA_APF, s_idx, s_rank] = deal(1,2,3,4,5,6,7,8,9,10,11,12);

data = csvread(slack_data, 1);

mpc = loadcase(casedata);
num_buses = size(mpc.bus, 1);
gen_len = size(mpc.gen, 1);


initial_res = problem(casedata);
for i = 1:num_buses
    if mpc.bus(i,BUS_TYPE) == 3
        slack_idx = i;
        slack_bus = mpc.bus(i, BUS_I);
        break;
    end
end

for idx = 1:gen_len
    if mpc.gen(idx, GEN_BUS) == slack_bus
        slack_gen = idx;
    end
end
slack_pg = initial_res.gen(slack_gen, PG);
fprintf(fileID, '%s, %s, %s, %s, %s, %s\n', 'BUS_I', 'BUS_TYPE', 'shed_bus', 'min_P_shed', 'num_conv', 'num_atmpt');
for b = 1:num_buses
    bus_num = data(b, s_BUS_I);
    mpc = loadcase(casedata);
    if mpc.bus(b, BUS_TYPE) ~= 1
        for i = 1:gen_len
            if mpc.gen(i,1) == bus_num
                gen_idx = i;
                break
            end
        end
    end 
    % set slackbus to PV, set slack power
    if b ~= slack_idx
        mpc.bus(slack_idx,BUS_TYPE) = 2;
        mpc.gen(slack_gen,PG) = slack_pg;
    end
    % add new gen for PQ buses
    if mpc.bus(b,BUS_TYPE) == 1
        mpc.gen(gen_len+1,:) = [mpc.bus(b,BUS_I),mpc.bus(b,PD),mpc.bus(b,QD),-1000,1000,0,100,1,1000,0,0,0,0,0,0,0,0,0,0,0,0];
        mpc.gencost(gen_len+1,:) = mpc.gencost(gen_len,:);
        gen_idx = gen_len + 1;
    end
    mpc.bus(b,BUS_TYPE) = 3;
    mpc.gen(gen_idx, VG) = data(b, s_VM);
    mpc.gen(gen_idx, VA) = 0;
    num_attempts = 0;
    num_conv = 0;
    shed_power_min = inf;
    for i = 1:num_buses
        if mpc.bus(i, PD) || mpc.bus(i, QD)
            num_attempts = num_attempts + 1;
            mpc_shed = mpc;
            mpc_shed.bus(i, [PD, QD]) = [0, 0];
            fprintf('Attack %i / %i: Shed %i / %i\n', b, num_buses, i, size(mpc_shed.bus, 1));
            results = runpf(mpc_shed, mpoption('verbose', 0, 'out.all', 0));
            if results.success == 1
                num_conv = num_conv + 1;
                shed_power = mpc.bus(i, PD);
                if abs(shed_power) < abs(shed_power_min)
                    shed_power_min = shed_power;
                    shed_bus = mpc_shed.bus(i, BUS_I);
                end
            end
        end
    end
    out = [data(b, 1), data(b, 2), shed_bus, shed_power_min, num_conv, num_attempts];
    fprintf(fileID, '%i, %i, %i, %i, %i, %i\n', out(1), out(2), out(3), out(4), out(5), out(6));
end
fclose(fileID);

