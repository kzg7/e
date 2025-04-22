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

