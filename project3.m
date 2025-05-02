% Active Suspension System: Pole Placement Control Design

mb = 300; mw = 60;
ks = 16000; kt = 190000;
bs = 1000;   

A = [  0,        1,         0,         0;
     -ks/mb, -bs/mb,  ks/mb,  bs/mb;
      0,        0,         0,         1;
      ks/mw,  bs/mw, -(ks + kt)/mw, -bs/mw ];

B = [0;
     1/mb;
     0;
    -1/mw];

C = [1 0 0 0 ; 0 0 1 0];     
D = 0; 
CAB = ctrb(A, B);  
rank(CAB)
OAC=obsv(A,C);
rank(OAC)
CAC_dual = ctrb(A.',C.');
rank(CAC_dual)
eig(A)
%Pole placement using Ackermann’s Formula
p1 = conv([1,5],[1,6]);
conv(p1,conv([1 8],[1 10]))
phiA = A^4 + 29*A^3 + 308*A^2 + 1420*A + 2400*eye(4);
KN = [0 0 0 1]*inv(CAB)*phiA
eig(A-B*KN)

desired_poles = [-5 -6 -8 -10];
K = place(A, B, desired_poles)
eig(A-B*K)

observer_poles = [-20 -25 -30 -35];  % Faster estimation
L = place(A.', C.', observer_poles)
Ke=L'

sys = A - B*K;
sys_op = ss(A, B, C, D);

sys_cl = ss(sys, B, C, D);

x0 = [0.1; 0; -0.05; 0]; 
t = 0:0.01:5;        

[y, t, x] = initial(sys_cl, x0, t);
[yo, to, xo] = initial(sys_op, x0, t);

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

x = x0;
xhat = [1; 0; 1; 0];  % initial observer estimate

X = zeros(4, length(t));
Xhat = zeros(4, length(t));
U = zeros(1, length(t));

for k = 1:length(t)
 
    y = C * x;
    u = -K * xhat;
    dx = A * x + B * u;
    dxhat = A * xhat + B * u + Ke * (y - C * xhat);
    
    x = x + dx * 0.01;
    xhat = xhat + dxhat * 0.01;
    
    X(:,k) = x;
    Xhat(:,k) = xhat;
    U(k) = u;
end

figure;
for i = 1:4
    subplot(2,2,i);
    plot(t, X(i,:), 'b', t, Xhat(i,:), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel(['State x_' num2str(i)]);
    legend('Actual', 'Estimated');
    title(['x_' num2str(i) ' vs xhat' num2str(i)]);
    grid on;
end


B_obsv(:,1) = B;
B_obsv(:,2) = Ke(:,1);
B_obsv(:,3) = Ke(:,2);
B_obsv

F = [0; 0; 0; kt/mw];
Be(:,1) = B;
Be(:,2) = F;
Be