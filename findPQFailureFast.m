function [S_high, angle, success, out, mpc] = findPQFailureFast(busIdx, mpc, tau, numAngles, Pmax)
define_constants;
if nargin < 5
    Pmax = 100 * max(mpc.bus(:, PD));
    if nargin < 4
        numAngles = 12;
        if nargin < 3
            tau = Pmax / 1e5;
        end
    end
end

if mpc.bus(busIdx, BUS_TYPE) ~= 1
    error('Bus must be PQ type');
end

angleStep = 2*pi/numAngles;
angles = zeros(1, numAngles);
for i = 1:numAngles
    angles(i) = (i-1) * angleStep;
end

P = 0;
Q = 0;
S = max(mpc.bus(:, PD));
Pinit = mpc.bus(busIdx, PD);
Qinit = mpc.bus(busIdx, QD);
Sinit = sqrt(Pinit^2 + Qinit^2);

results = runpf(mpc, mpoption('verbose', 0, 'out.all', 0));
failed = 0;
S_low = 0;
while (failed==0)&& (S <= Pmax)
    S = 2 * S;
    for idx = 1:numAngles
        P = S * cos(angles(idx));
        Q = S * sin(angles(idx));
        mpc.bus(busIdx, PD) = Pinit + P;
        mpc.bus(busIdx, QD) = Qinit + Q;
        results = runpf(mpc, mpoption('verbose', 0, 'out.all', 0));
        if results.success == 0
            failed = 1;
            angle = angles(idx);
            break
        else
            S_low = S;
        end
    end
end
success = failed;
S_high = S;
S = (S_low + S_high)/2;
while (S_high - S) > tau
    failed = 0;
    for idx = 1:numAngles
        P = S * cos(angles(idx));
        Q = S * sin(angles(idx));
        mpc.bus(busIdx, PD) = Pinit + P;
        mpc.bus(busIdx, QD) = Qinit + Q;
        results = runpf(mpc, mpoption('verbose', 0, 'out.all', 0));
        if results.success == 0
            failed = 1;
            angle = angles(idx);
            break
        end
    end
    if (failed == 1)
        S_high = S;
        S = (S + S_low)/2;
    else
        S_low = S;
        S = (S + S_high)/2;
    end
    
end
%success = failed;
out = [busIdx, P, Q, S_high, success];

    
    