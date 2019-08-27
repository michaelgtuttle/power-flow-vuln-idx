function [out_file] = slacktest_full(case_name,problem,problem_name)

if nargin < 3
    problem_name = 'pf';
    if nargin < 2
        problem = @(x) runpf(x, mpoption('verbose', 0, 'out.all', 0));
    end
end


    define_constants;

    % generate file name
    file_name = sprintf('slacktest-full-%s-%s.csv',case_name,problem_name);
    fileID = fopen(file_name,'w');
    
    mpc = loadcase(case_name);
    num_buses = size(mpc.bus,1);
    for i = 1:num_buses
        if mpc.bus(i,BUS_TYPE) == 3
            slack_idx = i;
            slack_bus = mpc.bus(i, BUS_I);
            break;
        end
    end
    
    
    initial_res = problem(case_name);
    num_gen = size(mpc.gen, 1);
    for idx = 1:num_gen
        if mpc.gen(idx, GEN_BUS) == slack_bus
            slack_gen = idx;
        end
    end
    
    slack_pg = initial_res.gen(slack_gen, PG);
    
    %bus, bus type, bus PD, bus QD, gen PG, gen QG, vm before failure, vm at failure, va at failure, va at failure
    fprintf(fileID,'BUS_I, BUS_TYPE, PD, QD, PG, QG, VM_BF, VM_AF, VA_BPF, VA_APF\n');
    
    for b = 1:num_buses
        fprintf('%i/%i\n\n', b, num_buses);
        go = 1;
        dec = 0.01;
        
        while go == 1
            mpc = loadcase(case_name);
            
            if b ~= slack_idx
                mpc.bus(slack_idx,BUS_TYPE) = 2;
                mpc.gen(slack_gen,PG) = slack_pg;
            end
            
            bus_num = mpc.bus(b, BUS_I);
            correctv = mpc.bus(b,VM);
            mpc.bus(b,VM) = mpc.bus(b,VM) - dec;
            faila = mpc.bus(b,VA);
            failv = mpc.bus(b,VM);
            gen_len = size(mpc.gen,1);
            bus_type_orig = mpc.bus(b,BUS_TYPE);
            
            
            if mpc.bus(b,BUS_TYPE) == 1
                mpc.gen(gen_len+1,:) = [mpc.bus(b,BUS_I),mpc.bus(b,PD),mpc.bus(b,QD),-1000,1000,failv,100,1,1000,0,0,0,0,0,0,0,0,0,0,0,0];
                mpc.gencost(gen_len+1,:) = mpc.gencost(gen_len,:);
                gen_idx = gen_len + 1;
                
            else
                for i = 1:gen_len
                    if mpc.gen(i,1) == bus_num
                        gen_idx = i;
                        mpc.gen(gen_idx,VG) = failv;
                        
                    end
                end
            end
             
            
            mpc.bus(b,BUS_TYPE) = 3;
            
            if failv <= 0
                % Does not consider voltages below 0.
                % if powerflow converges at 0 V, the algorithm is not
                % sensitive to this bus voltage
                go = 0;
                d_BUS_I = mpc.bus(b, BUS_I);
                d_PD = results.bus(b, PD);
                d_QD = results.bus(b, QD);
                disp(d_BUS_I)
                disp(bus_type_orig)
                disp(gen_idx)
                d_PG = results.gen(gen_idx, PG);
                d_QG = results.gen(gen_idx, QG);
                % append results to file
                % format: bus, bus type, P demand, Q demand, P generated, Q generated, vm before failure, 0 Volts (does not check V below 0), va at failure, va at failure
                fprintf(fileID,'%i, %i, %i, %i, %i, %i, %i, %i, %i, %i\n',d_BUS_I,bus_type_orig,d_PD,d_QD,d_PG,d_QG,correctv,0,faila,faila);
            
            else
                results = problem(mpc);
                dec = dec + 0.01;
            end
            
            
            if failv > 0 && results.success == 0
                go = 0;
                % append results to file
                % format: bus, bus type, bus PD, bus QD, gen PG, gen QG, vm before failure, vm at failure, va at failure, va at failure
                % after power flow
                fprintf(fileID,'%i, %i, %i, %i, %i, %i, %i, %i, %i, %i\n',mpc.bus(b,BUS_I),bus_type_orig,results.bus(b,PD),results.bus(b,QD),results.gen(gen_idx,PG),results.gen(gen_idx,QG),correctv,failv,results.bus(b,VA),faila);
                
            end
        end
    end
    
    fclose(fileID);
    out_file = file_name;
end
