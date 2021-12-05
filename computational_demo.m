% All values are in SI base units.
global STD_GRAVITY
global COULOMB_CONSTANT
global CHARGE1
global CHARGE2
global MASS

STD_GRAVITY = 9.80665;
COULOMB_CONSTANT = 8.9875517923e9;
CHARGE1 = 2.13e-07;
CHARGE2 = 2.17e-07;
MASS = 7.61e-04;
HEIGHT0 = 1.0;
TIMESTEP = 0.01;
NSTEP = 100;

% Initialize some variables
y = HEIGHT0;
vy = 0.0;
ys = 1: NSTEP;
vys = 1: NSTEP;

% Leapfrog integration
ay = force_y(y) / MASS;
for istep = 1: NSTEP
    % Leapfrog algorithm for a
    % single step.
    vy_half = vy + TIMESTEP * ay / 2;
    y_next = y + TIMESTEP * vy_half;
    ay_next = force_y(y_next) / MASS;
    vy_next = vy_half + TIMESTEP * ay_next / 2;

    % Store results.
    ys(istep) = y_next;
    vys(istep) = vy_next;

    % Move "next" results to current.
    y = y_next;
    vy = vy_next;
    ay = ay_next;
end

% Plot of position, velocity and energies as a function of time.

ts = (1 : NSTEP) * TIMESTEP;

hold on;
subplot(1, 3, 1)
plot(ts, ys)
xlabel("time t [s]")
ylabel("height y [m]")

subplot(1, 3, 2)
plot(ts, vys)
xlabel("time t [s]")
ylabel("velocity v_y [m/s]")

subplot(1, 3, 3)
eks = kinetic_energy(vys);
eps = potential_energy(ys);
plot(ts, eks, "b", ts, eps, "r", ts, eks + eps, "k")
%plot(ts, eps)
%plot(ts, eks + eps, 'color', "k")
legend("kin. E", "pot. E", "total E")
xlabel("time t [s]")
ylabel("energy [J]")
saveas(gca(), "matlab_vel_pos_ener.png", "png")
clf()
hold off;


% Plot velocity versus position

hold on;
plot(vys, ys)
xlabel("velocity v_y [m/s]")
ylabel("height y [m]")
xlim([-3.5 3.5])
ylim([0.02 1.05])
saveas(gca(), "matlab_vel_v_pos.png", "png")
clf()
hold off;


% Plot contour lines for the total energy

hold on;
[vy_grid, y_grid] = meshgrid(linspace(-3.5, 3.5, 51), linspace(0.02, 1, 51));
energy_grid = kinetic_energy(vy_grid) + potential_energy(y_grid);
contour(vy_grid, y_grid, energy_grid, 40)
xlabel("velocity v_y [m/s]")
ylabel("height y [m]")
saveas(gca(), "matlab_ener_contour.png", "png")
clf()
hold off;

function [out] = force_y(y)
    global MASS
    global STD_GRAVITY
    out = -MASS * STD_GRAVITY;
end

function [out] = potential_energy(y)
    global MASS
    global STD_GRAVITY
    out = MASS * STD_GRAVITY * y;
end

function [out] = kinetic_energy(vy)
    global MASS
    out = (0.5 * MASS) * vy.^2;
end
