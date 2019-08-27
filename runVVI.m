function [VVI, out] = runVVI(mpc, write)

if nargin < 2
    write = 0
end

define_constants;
slack = zeros(3, 1);
slack(1) = find(mpc.bus(:, BUS_TYPE) == 3);
slackBus = mpc.bus(slack(1), BUS_I);
slack(2) = find(mpc.gen(:, GEN_BUS) == slackBus);
results = runpf(mpc, mpoption('verbose', 0, 'out.all', 0));
slack(3) = results.gen(slack(2), PG);

numBuses = size(mpc.bus, 1);
VVI = zeros(numBuses, 1);
out = zeros(numBuses, 6);

for i = 1:numBuses
    [VVI(i), out(i, :), ~] = findVoltageFailure(i, mpc, slack);
end

if write == 1
    filename = sprintf('runVVI-%ibus.csv', numBuses);
    csvwrite(filename, out);
end
