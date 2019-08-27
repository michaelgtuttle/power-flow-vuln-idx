function [S, angles, PQBuses, success] = runPVI(mpc, tau, numAngles, Pmax, write)

define_constants;
if ischar(mpc)
    mpc = loadcase(mpc);
end
PQBuses = find(mpc.bus(:, BUS_TYPE) == 1);
if nargin < 5
    write = 1;
    if nargin < 4
        Pmax = 1000 * max(mpc.bus(:, PD));
        if nargin < 3
            numAngles = 12;
            if nargin < 2
                tau = max(mpc.bus(PQBuses, PD)) / 1000;
            end
        end
    end
end



PQBuses = find(mpc.bus(:, BUS_TYPE) == 1);
numPQBuses = size(PQBuses, 1);
PVI = zeros(numPQBuses, 2);
out = zeros(numPQBuses, 5);
success = zeros(numPQBuses, 1);
for idx = 1:numPQBuses
    bus = PQBuses(idx);
    [PVI(idx, 1), PVI(idx, 2), success(idx), out(idx, :)] = findPQFailureFast(bus, mpc, tau, numAngles, Pmax);
    if success == 0
        [PVI(idx, :), out(idx, :)] = deal(inf);
    end
end

S = PVI(:, 1);
angles = PVI(:, 2);
if write == 1
    fname = sprintf('runPVI-%ibus.csv', size(mpc.bus, 1));
    csvwrite(fname, out);
end
