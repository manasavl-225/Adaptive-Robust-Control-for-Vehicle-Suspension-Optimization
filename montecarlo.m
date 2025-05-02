
%% --- Nominal Parameters ---
mb_nom = 300; mw_nom = 60;
ks_nom = 16000; kt_nom = 190000;
bs_nom = 1000;

%% --- Simulation Setup ---
N = 100;  % Monte Carlo runs
t = 0:0.01:5;
x0 = [0.1; 0; -0.05; 0];  % Initial condition

% Preallocate storage
x_b_all = zeros(length(t), N);
x_w_all = zeros(length(t), N);
acc_b_all = zeros(length(t), N);
acc_w_all = zeros(length(t), N);

% Observer storage for first run only (for plotting)
x_true_hist = zeros(length(t), 4);
x_hat_hist = zeros(length(t), 4);

for i = 1:N
    % --- Parameter perturbation (±10%) ---
    mb = mb_nom * (1 + 0.2*(rand - 0.5));
    mw = mw_nom * (1 + 0.2*(rand - 0.5));
    ks = ks_nom * (1 + 0.2*(rand - 0.5));
    kt = kt_nom * (1 + 0.2*(rand - 0.5));
    bs = bs_nom * (1 + 0.2*(rand - 0.5));

    % --- State-space matrices ---
    A = [  0,        1,         0,         0;
         -ks/mb, -bs/mb,  ks/mb,  bs/mb;
          0,        0,         0,         1;
          ks/mw,  bs/mw, -(ks + kt)/mw, -bs/mw ];

    B = [0; 1/mb; 0; -1/mw];
    C = eye(4);
    D = zeros(4,1);

    % --- Controller Pole Placement ---
    desired_poles = [-5 -6 -8 -10];
    K = place(A, B, desired_poles);
    A_cl = A - B*K;

    % --- Observer Design ---
    observer_poles = [-20 -25 -30 -35];   % faster than controller
    Ke = place(A', C', observer_poles)';
    L = Ke';

    % --- Initialize states ---
    x_true = x0;         % Actual system state
    x_hat = zeros(4,1);  % Estimated state

    x_b = zeros(length(t),1);
    x_w = zeros(length(t),1);
    acc_b = zeros(length(t),1);
    acc_w = zeros(length(t),1);

    for j = 1:length(t)
        % Control input based on estimated state
        u = -K * x_hat;

        % System dynamics (Euler integration)
        dx = A * x_true + B * u;
        x_true = x_true + dx * 0.01;

        % Observer dynamics
        y = C * x_true;
        y_hat = C * x_hat;
        dx_hat = A * x_hat + B * u + L * (y - y_hat);
        x_hat = x_hat + dx_hat * 0.01;

        % Store body & wheel displacements
        x_b(j) = x_true(1);
        x_w(j) = x_true(3);

        % Compute accelerations
        acc_b(j) = (-ks*(x_true(1) - x_true(3)) - bs*(x_true(2) - x_true(4))) / mb;
        acc_w(j) = ( ks*(x_true(1) - x_true(3)) + bs*(x_true(2) - x_true(4)) - kt*x_true(3) ) / mw;

        % Save for plotting if first run
        if i == 1
            x_true_hist(j,:) = x_true';
            x_hat_hist(j,:) = x_hat';
        end
    end

    % Store for Monte Carlo results
    x_b_all(:,i) = x_b;
    x_w_all(:,i) = x_w;
    acc_b_all(:,i) = acc_b;
    acc_w_all(:,i) = acc_w;
end

%% --- Monte Carlo Plotting: Displacement & Acceleration ---
figure;
subplot(2,2,1);
plot(t, x_b_all, 'Color', [0.5 0.5 1]); hold on;
plot(t, mean(x_b_all,2), 'b', 'LineWidth', 2);
title('Monte Carlo: x_b (m)'); ylabel('x_b (m)'); xlabel('Time (s)'); grid on;

subplot(2,2,2);
plot(t, x_w_all, 'Color', [0.5 0.5 1]); hold on;
plot(t, mean(x_w_all,2), 'b', 'LineWidth', 2);
title('Monte Carlo: x_w (m)'); ylabel('x_w (m)'); xlabel('Time (s)'); grid on;

subplot(2,2,3);
plot(t, acc_b_all, 'Color', [1 0.6 0.6]); hold on;
plot(t, mean(acc_b_all,2), 'r', 'LineWidth', 2);
title('Monte Carlo: \ddot{x}_b (m/s^2)'); ylabel('\ddot{x}_b'); xlabel('Time (s)'); grid on;

subplot(2,2,4);
plot(t, acc_w_all, 'Color', [1 0.6 0.6]); hold on;
plot(t, mean(acc_w_all,2), 'r', 'LineWidth', 2);
title('Monte Carlo: \ddot{x}_w (m/s^2)'); ylabel('\ddot{x}_w'); xlabel('Time (s)'); grid on;

%% --- Observer Performance (First Simulation) ---
figure;
for k = 1:4
    subplot(2,2,k);
    plot(t, x_true_hist(:,k), 'k', 'LineWidth', 1.5); hold on;
    plot(t, x_hat_hist(:,k), '--r', 'LineWidth', 1.5);
    legend('True', 'Estimated');
    title(['State x_', num2str(k), ' vs \hat{x}_', num2str(k)]);
    xlabel('Time (s)');
    grid on;
end

figure;
for k = 1:4
    subplot(2,2,k);
    plot(t, x_true_hist(:,k) - x_hat_hist(:,k), 'b', 'LineWidth', 1.5);
    title(['Estimation Error: x_', num2str(k), ' - \hat{x}_', num2str(k)]);
    xlabel('Time (s)');
    ylabel('Error');
    grid on;
end
