function [MVAbase, bus, gen, branch, success, et, converged, Vlf, Sbuslf] = my_fdi(casedata, injections, mpopt, fname, solvedcase)
% injections = {{bus_num1, P + jQ, 'S'}, {bus_num2, V * exp(j*theta *
% pi/180), 'V'}, ...}


%RUNSE  Runs a state estimator.
%   [BASEMVA, BUS, GEN, BRANCH, SUCCESS, ET] = ...
%           RUNSE(CASEDATA, MPOPT, FNAME, SOLVEDCASE)
%
%   Runs a state estimator (after a Newton power flow). Under construction with
%   parts based on code from James S. Thorp.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   parts based on code by James S. Thorp, June 2004
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialize  -----
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% default arguments
if nargin < 5
    solvedcase = '';                %% don't save solved case
    if nargin < 4
        fname = '';                 %% don't print results to a file
        if nargin < 3
            mpopt = mpoption('verbose', 0, 'out.all', 0); 
            %% use default options
            if nargin < 2
                injections = {};
                if nargin < 1
                    casedata = 'case9'; %% default data file is 'case9.m'
                end
            end
        end
    end
end

%% options
qlim = mpopt.pf.enforce_q_lims;         %% enforce Q limits on gens?
dc = strcmp(upper(mpopt.model), 'DC');  %% use DC formulation?

%% read data
mpc = loadcase(casedata);

%% add zero columns to branch for flows if needed
if size(mpc.branch,2) < QT
  mpc.branch = [ mpc.branch zeros(size(mpc.branch, 1), QT-size(mpc.branch,2)) ];
end

%% convert to internal indexing
[baseMVA, bus, gen, branch] = loadcase(casedata);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
mpc = ext2int(mpc);
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

%%-----  run the power flow  -----
t0 = clock;
its = 0;            %% total iterations
if mpopt.verbose > 0
    v = mpver('all');
    fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
end
if dc                               %% DC formulation
    if mpopt.verbose > 0
      fprintf(' -- DC Power Flow\n');
    end
    %% initial state
    Va0 = bus(:, VA) * (pi/180);
    
    %% build B matrices and phase shift injections
    [B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
    
    %% compute complex bus power injections (generation - load)
    %% adjusted for phase shifters and real shunts
    Pbus = real(makeSbus(baseMVA, bus, gen)) - Pbusinj - bus(:, GS) / baseMVA;
    
    %% "run" the power flow
    [Va, F] = dcpf(B, Pbus, Va0, ref, pv, pq);
    its = 1;
    
    %% update data matrices with solution
    branch(:, [QF, QT]) = zeros(size(branch, 1), 2);
    branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
    branch(:, PT) = -branch(:, PF);
    bus(:, VM) = ones(size(bus, 1), 1);
    bus(:, VA) = Va * (180/pi);
    %% update Pg for slack generator (1st gen at ref bus)
    %% (note: other gens at ref bus are accounted for in Pbus)
    %%      Pg = Pinj + Pload + Gs
    %%      newPg = oldPg + newPinj - oldPinj
    refgen = zeros(size(ref));
    for k = 1:length(ref)
        temp = find(gbus == ref(k));
        refgen(k) = on(temp(1));
    end
    gen(refgen, PG) = gen(refgen, PG) + (B(ref, :) * Va - Pbus(ref)) * baseMVA;
else                                %% AC formulation
    alg = upper(mpopt.pf.alg);
    if mpopt.verbose > 0
        switch alg
            case 'NR'
                solver = 'Newton';
            case 'FDXB'
                solver = 'fast-decoupled, XB';
            case 'FDBX'
                solver = 'fast-decoupled, BX';
            case 'GS'
                solver = 'Gauss-Seidel';
            otherwise
                solver = 'unknown';
        end
        fprintf(' -- AC Power Flow (%s)\n', solver);
    end
    %% initial state
    % V0    = ones(size(bus, 1), 1);            %% flat start
    V0  = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
    vcb = ones(size(V0));           %% create mask of voltage-controlled buses
    vcb(pq) = 0;                    %% exclude PQ buses
    k = find(vcb(gbus));            %% in-service gens at v-c buses
    V0(gbus(k)) = gen(on(k), VG) ./ abs(V0(gbus(k))).* V0(gbus(k));
    
    if qlim
        ref0 = ref;                         %% save index and angle of
        Varef0 = bus(ref0, VA);             %%   original reference bus(es)
        limited = [];                       %% list of indices of gens @ Q lims
        fixedQg = zeros(size(gen, 1), 1);   %% Qg of gens at Q limits
    end

    %% build admittance matrices
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
    
    repeat = 1;
    while (repeat)
        %% function for computing V dependent complex bus power injections
        %% (generation - load)
        Sbus = @(Vm)makeSbus(baseMVA, bus, gen, mpopt, Vm);
        
        %% run the power flow
        switch alg
            case 'NR'
                [V, success, iterations] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
            case {'FDXB', 'FDBX'}
                [Bp, Bpp] = makeB(baseMVA, bus, branch, alg);
                [V, success, iterations] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt);
            case 'GS'
                if (~isempty(mpopt.exp.sys_wide_zip_loads.pw) && ...
                        any(mpopt.exp.sys_wide_zip_loads.pw(2:3))) || ...
                        (~isempty(mpopt.exp.sys_wide_zip_loads.qw) && ...
                        any(mpopt.exp.sys_wide_zip_loads.qw(2:3)))
                    warning('runpf: Gauss-Seidel algorithm does not support ZIP load model. Converting to constant power loads.')
                    mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', ...
                                    struct('pw', [], 'qw', []));
                end
                [V, success, iterations] = gausspf(Ybus, Sbus([]), V0, ref, pv, pq, mpopt);
            otherwise
                error('runpf: Only Newton''s method, fast-decoupled, and Gauss-Seidel power flow algorithms currently implemented.');
        end
        its = its + iterations;
        repeat = 0;
    end    
        %% update data matrices with solution
        [bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
end
et = etime(clock, t0);

%%--------------------  begin state estimator code  --------------------
%% save some values from load flow solution
Pflf=branch(:,PF);
Qflf=branch(:,QF);
Ptlf=branch(:,PT);
Qtlf=branch(:,QT);
Sbuslf = V .* conj(Ybus * V);
Vlf=V;
for i = 1:size(injections, 2)
    attack = injections{i};
    bus_index = find(bus(:, BUS_I) == attack{1});
    attack_val = attack{2};
    if attack{3} == 'S'
        Sbuslf(bus_index) = Sbuslf(bus_index) + attack_val / baseMVA;
    elseif attack{3} == 'V'
        Vlf(bus_index) = attack_val;
    end
end

    

%% run state estimator
[V, converged, i] = state_est(branch, Ybus, Yf, Yt, Sbuslf, Vlf, ref, pv, pq, mpopt);

%% update data matrices to match estimator solution ...
%% ... bus injections at PQ buses
Sbus = V .* conj(Ybus * V);
bus(pq, PD) = -real(Sbus(pq)) * baseMVA;
bus(pq, QD) = -imag(Sbus(pq)) * baseMVA;
%% ... gen outputs at PV buses
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?
gen(on, PG) = real(Sbus(gbus)) * baseMVA + bus(gbus, PD);   %% inj P + local Pd
%% ... line flows, reference bus injections, etc.
[bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);

%% plot differences from load flow solution
Pfe=branch(:,PF);
Qfe=branch(:,QF);
Pte=branch(:,PT);
Qte=branch(:,QT);
nbr = length(Pfe);
subplot(3,2,1), plot(180/pi*(angle(Vlf)-angle(V)),'.'), title('Voltage Angle (deg)');
subplot(3,2,2), plot(abs(Vlf)-abs(V),'.'), title('Voltage Magnitude (p.u.)');
subplot(3,2,3), plot((1:nbr),(Pfe-Pflf),'r.',(1:nbr),(Pte-Ptlf),'b.'), title('Real Flow (MW)');
subplot(3,2,4), plot((1:nbr),(Qfe-Qflf),'r.',(1:nbr),(Qte-Qtlf),'b.'), title('Reactive Flow (MVAr)');
subplot(3,2,5), plot(baseMVA*real(Sbuslf-Sbus), '.'), title('Real Injection (MW)');
subplot(3,2,6), plot(baseMVA*imag(Sbuslf-Sbus), '.'), title('Reactive Injection (MVAr)');
%%--------------------  end state estimator code  --------------------

%%-----  output results  -----
%% convert back to original bus numbering & print results
[bus, gen, branch] = int2ext(i2e, bus, gen, branch);
if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        if mpopt.out.all == 0
            printpf(baseMVA, bus, gen, branch, [], success, et, fd, ...
                mpoption(mpopt, 'out.all', -1));
        else
            printpf(baseMVA, bus, gen, branch, [], success, et, fd, mpopt);
        end
        fclose(fd);
    end
end
printpf(baseMVA, bus, gen, branch, [], success, et, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, baseMVA, bus, gen, branch);
end

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end
