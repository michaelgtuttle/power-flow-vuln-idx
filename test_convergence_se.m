function [] = test_convergence_se(casedata, radial_step)
define_constants;
mpc = loadcase(casedata);
num_buses = size(mpc.bus, 1);
Smag = 0;
theta = 0;
theta_inc = pi / 6;
ang = 0;
filename = sprintf('test-convergence-se-%s-V.csv', casedata);
filename2 = sprintf('test-convergence-se-%s-S.csv', casedata);
fileID = fopen(filename, 'w');
fileID2 = fopen(filename2, 'w');
[~, ~, ~, ~, ~, ~, ~, Vlf, Sbuslf] = my_fdi('case30', {});
Sbus_mag = abs(Sbuslf);
Sbus_theta = atan(imag(Sbuslf) ./ real(Sbuslf));
V_mag = abs(Vlf);
V_theta = atan(imag(Vlf) ./ real(Vlf));
for b = 1:num_buses
    bus_num = mpc.bus(b, BUS_I);
    Smag = 0;
    theta = 0;
    failed = 0;
    V = V_mag(b);
    ang = 0;
    %[~,~,~,~,~,~, converged] = my_fdi(casedata);
    while 1
        fprintf('%i/%i: %i < %.1f\n', b, num_buses, V, ang);
        if mpc.bus(b, BUS_TYPE) == 3
            success = 0;
            fprintf(fileID, '%i, %i, %i, %i, %i, %i\n', bus_num, V, ang + V_theta(b)*180/pi, V_mag(b), V_theta(b), success);
            break
        end
        injections = {{bus_num, V*exp(1j*(V_theta(b) + ang*pi/180)), 'V'}};
        [~,~,~,~,~,~, converged] = my_fdi(casedata, injections);
        if ~converged
            success = 1;
            fprintf(fileID, '%i, %i, %i, %i, %i, %i\n', bus_num, V, ang + V_theta(b)*180/pi, V_mag(b), V_theta(b), success);
            break
        end
        if V < 0
            success = 0;
            fprintf(fileID, '%i, %i, %i, %i, %i, %i\n', bus_num, V, ang + V_theta(b)*180/pi, V_mag(b), V_theta(b), success);
            break
        end
        ang = -1 * ang;
        if ang >= 0
            ang = ang + 2.5;
        end
        if ang > 10
            ang = 0;
            V = V - 0.01;
        end
    end
end

        
for b = 1:num_buses
    bus_num = mpc.bus(b, BUS_I);
    Smag = 0;
    theta = 0;
    failed = 0;
    V = V_mag(b);
    ang = 0;
    while 1
        if mpc.bus(b, BUS_TYPE) == 3
            success = 0;
            fprintf(fileID2, '%i, %i, %i, %i, %i, %i\n', bus_num, Smag, theta, success, Sbus_mag(b), Sbus_theta(b));
            break
        end
        fprintf('BUS %i/%i: S = %i <%i\n', b, num_buses, Smag, theta * 180/pi);
        S = Smag * exp(1j * theta);
        injections = {{bus_num, S, 'S'}};
        [~,~,~,~,~,~, converged] = my_fdi(casedata, injections);
        if ~converged
            success = 1;
            failed = 1;
            fprintf(fileID2, '%i, %i, %i, %i, %i, %i\n', bus_num, Smag, theta, success, Sbus_mag(b), Sbus_theta(b));
            break;
        end
        theta = theta + theta_inc;
        if abs(theta - 2*pi) < 0.001
            if failed == 1
                break;
            end
            theta = 0;
            Smag = Smag + radial_step;
        end
    end
end
fclose(fileID);
fclose(fileID2);
end