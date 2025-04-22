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
