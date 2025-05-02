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

F = [0; 0; 0 ;kt/mw]; 
CAB = ctrb(A, B); 
rank(CAB)
OAC=obsv(A,C);
rank(OAC)
CAC_dual = ctrb(A.',C.');
rank(CAC_dual)


Q = diag([100, 1, 100, 1]);  
R = 1;                       
Kqr = lqr(A, B, Q, R)
[Mbar,Kbar,L] = icare(A,B,Q,R)
eig(A-B*Kbar)

F = [0; 0; 0; kt/mw];    
V = 0.01;                
W = 0.001 * eye(2);      

Qk = F * V * F.';       

[P, L_raw, ~] = icare(A', C', Qk, W);

Ke = L_raw';  

disp('Kalman gain :');
disp(Ke);

T = 0:0.001:5;
x = [0.05; 0; -0.02; 0];       
xhat = zeros(4,1);             

w = 0.01 * sin(2*pi*1*T);      % Road disturbance
noise = 0.001*randn(2, length(T));  % Sensor noise

X = zeros(4, length(T));
Xhat = zeros(4, length(T));
acc = zeros(1,length(T));
U = zeros(1,length(T));

dt = T(2) - T(1);

for k = 1:length(T)
    y = C*x + noise(:,k);             
    yhat = C*xhat;

    u = -Kqr * xhat;                  
    dx = A*x + B*u + F*w(k);         
    dxhat = A*xhat + B*u + Ke*(y - yhat);  

    x = x + dx*dt;
    xhat = xhat + dxhat*dt;

    X(:,k) = x;
    Xhat(:,k) = xhat;
    acc(k) = dx(2);
    U(k) = u;
end


figure;

subplot(4,1,1)
plot(T, X(1,:), 'b', 'LineWidth', 1.5)
title('Body Displacement (m)'); ylabel('x_b (m)'); grid on;

subplot(4,1,2)
plot(T, X(3,:), 'r', 'LineWidth', 1.5)
title('Wheel Displacement (m)'); ylabel('x_w (m)'); grid on;

subplot(4,1,3)
plot(T, acc, 'g', 'LineWidth', 1.5)
title('Body Acceleration (m/s^2)'); ylabel('a_b'); grid on;

subplot(4,1,4)
plot(T, w, 'k--', 'LineWidth', 1.5)
title('Road Disturbance Input'); xlabel('Time (s)'); ylabel('w (m)'); grid on;

sgtitle('Active Suspension with LQR + Kalman Filter');


B_obsv(:,1) = B;
B_obsv(:,2) = Ke(:,1);
B_obsv(:,3) = Ke(:,2);
B_obsv