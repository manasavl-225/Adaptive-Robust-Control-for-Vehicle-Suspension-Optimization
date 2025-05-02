% Active Suspension System: Pole Placement Control Design
%% --- Parameters ---
mb = 300; mw = 60;
ks = 16000; kt = 190000;
bs = 1000;   % damping coefficient (N·s/m)

%% --- State-Space Matrices ---
% State vector: x = [x_b; ẋ_b; x_w; ẋ_w]
A = [  0,        1,         0,         0;
     -ks/mb, -bs/mb,  ks/mb,  bs/mb;
      0,        0,         0,         1;
      ks/mw,  bs/mw, -(ks + kt)/mw, -bs/mw ];

B = [0;
     1/mb;
     0;
    -1/mw];

C = [1 0 0 0 ; 0 0 1 0];     % Full state output
D = 0; % No direct feedthrough
CAB = ctrb(A, B);  % use only the actuator input (f_s)
rank(CAB)
OAC=obsv(A,C);
rank(OAC)
CAC_dual = ctrb(A.',C.');
rank(CAC_dual)
eig(A)
p1 = conv([1,5],[1,6]);
conv(p1,conv([1 8],[1 10]))
phiA = A^4 + 29*A^3 + 308*A^2 + 1420*A + 2400*eye(4);
KN = [0 0 0 1]*inv(CAB)*phiA
eig(A-B*KN)

desired_poles = [-5 -6 -8 -10];
K = place(A, B, desired_poles) %check
eig(A-B*K)

%designing the observer
observer_poles = [-20 -25 -30 -35];  % Faster estimation
L = place(A.', C.', observer_poles)
Ke=L';

% Closed-loop system (A - B*K)
sys = A - B*K;
sys_op = ss(A, B, C, D);
% Create closed-loop state-space model
sys_cl = ss(sys, B, C-D*K, D);

%% --- Simulation ---
% Initial conditions: e.g., body displacement 0.1 m, wheel -0.05 m
x0 = [0.1; 0; -0.05; 0];  % [x_b; ẋ_b; x_w; ẋ_w]
t = 0:0.01:5;             % simulation time

% Simulate response to initial condition
[y, t, x] = initial(sys_cl, x0, t);
[yo, to, xo] = initial(sys_op, x0, t);

%% --- Plot Results ---
figure;
plot(t, x, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('State values');
legend('x_b (m)', 'ẋ_b (m/s)', 'x_w (m)', 'ẋ_w (m/s)');
title('with feedback State Response with Pole Placement');
grid on;
figure;
plot(to, xo, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('State values');
legend('x_b (m)', 'ẋ_b (m/s)', 'x_w (m)', 'ẋ_w (m/s)');
title('without feedback State Response');
grid on;


B_obsv(:,1) = B;
B_obsv(:,2) = Ke(:,1);
B_obsv(:,3) = Ke(:,2);
B_obsv

