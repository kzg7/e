% mag_field.m
% M-file to calculate and plot the net magnetic field produced by a three-phase stator.

% Set up basic conditions
bmax = 1;               % Normalize maximum magnetic field to 1
freq = 60;              % Frequency in Hz
w = 2 * pi * freq;      % Angular velocity (rad/s)

% Time vector for one cycle (1/60 s) with 6000 points
t = 0:1/6000:1/60;

% Generate the three component magnetic fields (complex representation)
% Baa, Bbb, Bcc are phasors with 120-degree phase differences
Baa = bmax * sin(w * t) .* (cos(0) + 1i * sin(0));                     % Phase A
Bbb = bmax * sin(w * t - 2 * pi / 3) .* (cos(2 * pi / 3) + 1i * sin(2 * pi / 3)); % Phase B
Bcc = bmax * sin(w * t + 2 * pi / 3) .* (cos(-2 * pi / 3) + 1i * sin(-2 * pi / 3)); % Phase C

% Calculate the net magnetic field
Bnet = Baa + Bbb + Bcc;

% Calculate a circle representing the expected maximum value of Bnet (radius = 1.5)
circle = 1.5 * (cos(w * t) + 1i * sin(w * t));

% Plot the magnitude and direction of the resulting magnetic fields
% Baa: black, Bbb: blue, Bcc: magenta, Bnet: red
figure;
for ii = 1:length(t)
    % Clear previous plot
    clf;
    
    % Plot the reference circle
    plot(circle, 'k--', 'LineWidth', 1); % Dashed black circle
    hold on;
    
    % Plot the four magnetic fields
    plot([0 real(Baa(ii))], [0 imag(Baa(ii))], 'k', 'LineWidth', 2); % Phase A
    plot([0 real(Bbb(ii))], [0 imag(Bbb(ii))], 'b', 'LineWidth', 2); % Phase B
    plot([0 real(Bcc(ii))], [0 imag(Bcc(ii))], 'm', 'LineWidth', 2); % Phase C
    plot([0 real(Bnet(ii))], [0 imag(Bnet(ii))], 'r', 'LineWidth', 3); % Net field
    
    % Set plot properties
    axis square;
    axis([-2 2 -2 2]);
    title('Magnetic Field Vectors in a Three-Phase Stator');
    xlabel('Real Part');
    ylabel('Imaginary Part');
    grid on;
    
    % Update plot
    drawnow;
    pause(0.01); % Small pause for animation effect
end
hold off;
