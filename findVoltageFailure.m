function [Vlow, out, failure] = findVoltageFailure(busIdx, mpc, slack)

define_constants;
gen_len = size(mpc.gen, 1);

vInitial = mpc.bus(busIdx, VM);


slackRow = mpc.gen(slack(2), :);
if busIdx ~= slack(1)
    mpc.bus(slack(1),BUS_TYPE) = 2;
    mpc.gen(slack(2),PG) = slack(3); %slack_gen power = slack_pg
end
busNum = mpc.bus(busIdx, BUS_I);

if mpc.bus(busIdx,BUS_TYPE) == 1
    mpc.gen(gen_len+1,:) = slackRow;
    %mpc.gen(gen_len+1,:) = [mpc.bus(busIdx,BUS_I),mpc.bus(busIdx,PD),mpc.bus(busIdx,QD),-1000,1000,1,100,1,1000,0,0,0,0,0,0,0,0,0,0,0,0];
    if(isfield(mpc, 'gencost'))
        mpc.gencost(gen_len+1,:) = mpc.gencost(gen_len,:);
    end
    genIdx = gen_len + 1;
    mpc.gen(genIdx, GEN_BUS) = busNum;
else    
    for i = 1:gen_len
        if mpc.gen(i,1) == busNum
            genIdx = i;
        end
    end
end

Vhigh = vInitial;
Vlow = 0;
V = 0.5;

mpc.bus(busIdx, BUS_TYPE) = 3;

failure = runpf(mpc, mpoption('verbose', 0, 'out.all', 0));

while (Vhigh - Vlow) >= 0.005 
    mpc.gen(genIdx, VG) = V;
    results = runpf(mpc, mpoption('verbose', 0, 'out.all', 0));
    if results.success == 0
        failure = results;
        Vlow = V;
        V = (V + Vhigh) / 2;
    else
        Vhigh = V;
        V = (V + Vlow) / 2;
    end
end

out = [busIdx, failure.bus(busIdx, PD), failure.gen(genIdx, PG), failure.bus(busIdx, QD), failure.gen(genIdx, QG), Vlow];
