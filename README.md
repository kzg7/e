%6A
% M-file, mag_ fi e l d.m
% M-file t o ca l c ulate the net magnetic fi e l d produ ced
% by a three-pha se s tat o r.
% Set up the basi c con d itio n s
bmax = 1; % No rmalize bmax t o 1
freq = 120; % 60 Hz
w = 2* pi * freq; % a n g ular ve l oc ity (rad/ s )
% Firs t , generate the three component magnetic fi e l ds
t = 0:1/6000:1/60;
Baa = sin(w*t).*(cos(0) + j*sin(0));
Bbb = sin(w*t - 2*pi/3) .* (cos (2*pi/3) + j*sin(2*pi/3)) ;
Bcc = sin(w*t + 2*pi/3) .* (cos (- 2*pi /3) + j*sin(- 2*pi /3)) ;

%Ca l c ulate Bnet
Bnet = Baa + Bbb + Bcc;

% Ca l c ulate a c irc l e representing the expected maximum
% va lue o f Enet

circle = 1.5 * (cos(w*t) + j *sin(w*t)) ;

% Plo t the magnitude and d irection of the r esulting magne ti c
% fi e l ds. No te that Baa i s b l ack, Bbb i s b lue, Bcc i s
% magenta , and Enet i s red.

for ii = 1:length(t)
    %   Plo t the reference c irc l e
    plot(circle,'k') ;
    hold on;
    % Plo t the f o ur magneti c fi e l ds
    plot ([0 real(Baa(ii))] , [ 0 imag(Baa(ii))] , 'k', 'LineWidth' ,2);
    plot ([0 real(Bbb(ii))] , [ 0 imag(Bbb(ii))] , 'b' , 'LineWidth' ,2) ;
    plot([0 real(Bcc(ii))] , [ 0 imag(Bcc(ii))] , 'm' , 'LineWidth' ,2) ;
    plot ( [0 real(Bnet(ii))] , [ 0 imag(Bnet (ii)) ] ,'r' ,'LineWidth',3 ) ;
    axis square;
    axis( [ -2 2 -2 2 ] ) ;
    drawnow;
    hold off ;
end


%6B
clc;
clear;
close all;

% Set up the basic conditions
bmax = 1; % Normalize bmax to 1
freq = 60; % 60 Hz
P = 2; % Number of poles (set to 2 for standard operation)
w = 2 * pi * freq * P / 2; % Adjusted angular velocity (pole dependent)

% Rotor parameters
rotor_speed = w / P; % Mechanical speed of rotor magnetic field
initial_angle = pi/4; % Initial misalignment angle

t = 0:1/6000:1/60; % Time vector

% Generate the three-phase magnetic fields
Baa = sin(w*t).*(cos(0) + j*sin(0));
Bbb = sin(w*t - 2*pi/3) .* (cos (2*pi/3) + j*sin(2*pi/3)) ;
Bcc = sin(w*t + 2*pi/3) .* (cos (- 2*pi /3) + j*sin(- 2*pi /3)) ;

% Calculate Bnet
Bnet = Baa + Bbb + Bcc;

% Rotor magnetic field (one-loop rotor interaction)
Brotor = bmax * (cos(rotor_speed * t + initial_angle) + j * sin(rotor_speed * t + initial_angle));

% Calculate a reference circle representing the expected maximum
circle = 1.5 * (cos(w*t) + j*sin(w*t));

% Plot the magnitude and direction of the resulting magnetic fields
figure;
for ii = 1:length(t)
    clf;
    % Plot the reference circle
    plot(circle, 'k');
    hold on;
    
    % Plot the three-phase stator magnetic fields
    plot([0 real(Baa(ii))], [0 imag(Baa(ii))], 'k', 'LineWidth', 2);
    plot([0 real(Bbb(ii))], [0 imag(Bbb(ii))], 'b', 'LineWidth', 2);
    plot([0 real(Bcc(ii))], [0 imag(Bcc(ii))], 'm', 'LineWidth', 2);
    
    % Plot the net rotating magnetic field
    plot([0 real(Bnet(ii))], [0 imag(Bnet(ii))], 'r', 'LineWidth', 3);
    
    % Plot the rotor magnetic field
    plot([0 real(Brotor(ii))], [0 imag(Brotor(ii))], 'g', 'LineWidth', 3);
    
    axis square;
    axis([-2 2 -2 2]);
    drawnow;
    hold off;
end


%7A

i_a = linspace(0, 60, 21);
e_a = 277.0;
x_s = 1.0;
pf_values = [0.2, 0.4, 0.6, 0.8];
v_phase_all = zeros(length(pf_values), length(i_a));

for p = 1:length(pf_values)
    pf = pf_values(p);
    theta = acos(pf);
    for ii = 1:length(i_a)
        I = i_a(ii);
        v_real = e_a - x_s * I * sin(theta);
        v_imag = x_s * I * cos(theta);
        v_phase_all(p, ii) = sqrt(v_real^2 + v_imag^2);
    end
end

v_t_all = v_phase_all * sqrt(3);

figure;
hold on;
colors = ['r', 'g', 'b', 'k'];
for p = 1:length(pf_values)
    plot(i_a, v_t_all(p, :), 'Color', colors(p), 'LineWidth', 2);
end
grid on;
xlabel('Armature Current (A)');
ylabel('Terminal Voltage (V)');
title('Terminal Characteristics of Synchronous Generator');
legend('PF = 0.2 lagging', 'PF = 0.4 lagging', 'PF = 0.6 lagging', 'PF = 0.8 lagging');


%7B
i_a = linspace(0, 60, 21);
e_a = 277.0;
x_s = 1.0;
pf_values = [0.2, 0.4, 0.6, 0.8];
v_phase_all = zeros(length(pf_values), length(i_a));

for p = 1:length(pf_values)
    pf = pf_values(p);
    theta = -acos(pf);
    for ii = 1:length(i_a)
        I = i_a(ii);
        v_real = e_a - x_s * I * sin(theta);
        v_imag = x_s * I * cos(theta);
        v_phase_all(p, ii) = sqrt(v_real^2 + v_imag^2);
    end
end

v_t_all = v_phase_all * sqrt(3);

figure;
hold on;
colors = ['r', 'g', 'b', 'k'];
for p = 1:length(pf_values)
    plot(i_a, v_t_all(p, :), 'Color', colors(p), 'LineWidth', 2);
end
grid on;
xlabel('Armature Current (A)');
ylabel('Terminal Voltage (V)');
title('Terminal Characteristics of Synchronous Generator (Leading PF)');
legend('PF = 0.2 leading', 'PF = 0.4 leading', 'PF = 0.6 leading', 'PF = 0.8 leading');

%8

r1 = 0.641;
x1 = 1.106;
r2 = 0.332;
x2 = 0.464;
xm = 26.3;
v_phase = 460 / sqrt(3);
n_sync = 1800;
w_sync = 188.5;

z_th = 1j * xm * (r1 + 1j * x1) / (r1 + 1j * (x1 + xm));
v_th = v_phase * (1j * xm) / (r1 + 1j * (x1 + xm));
r_th = real(z_th);
x_th = imag(z_th);

s = linspace(0.001, 1, 51);
n_m = (1 - s) * n_sync;

t_ind1 = zeros(1, 51);
t_ind2 = zeros(1, 51);
t_ind3 = zeros(1, 51);

for ii = 1:51
    t_ind1(ii) = (3 * abs(v_th)^2 * r2 / s(ii)) / ...
        (w_sync * ((r_th + r2 / s(ii))^2 + (x_th + x2)^2));
    t_ind2(ii) = (3 * abs(v_th)^2 * (2 * r2) / s(ii)) / ...
        (w_sync * ((r_th + (2 * r2) / s(ii))^2 + (x_th + x2)^2));
    t_ind3(ii) = (3 * abs(v_th)^2 * (0.5 * r2) / s(ii)) / ...
        (w_sync * ((r_th + (0.5 * r2) / s(ii))^2 + (x_th + x2)^2));
end

figure;
plot(n_m, t_ind1, 'k', 'LineWidth', 2.0);
hold on;
plot(n_m, t_ind2, 'k--', 'LineWidth', 2.0);
plot(n_m, t_ind3, 'k:', 'LineWidth', 2.0);
xlabel('n_m (rpm)', 'FontWeight', 'bold');
ylabel('\tau_{ind} (Nm)', 'FontWeight', 'bold');
title('Induction Motor Torque-Speed Characteristic', 'FontWeight', 'bold');
legend('Original R_2', '2 × R_2', '0.5 × R_2');
grid on;
hold off;

