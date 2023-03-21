clc; clear; close all;
tic

%% Syms
syms y [1 2] real
syms r [1 14] real
syms rr [1 2] real
syms rrr [1 12] real
y = y';
r1 = r(1 : 14)';
rr = rr';
rrr = rrr';

%% System & Parameter
alpha = 1.5;
d = 2;
Theta = [1, 0; 0, 1];

A{1} = [1 - alpha^2, -1; 0.45,  -0.1];
A{2} = [          1, -1; 0.45,  -0.1];

D{1} = [1; 0];  % v = 1
D{2} = [0; 2];  % v = 2
% B{1} = [1; 0];   % v = 1
% B{2} = [2; 0];   % v = 2

maxT1 = 1;
minT3 = 0;
I = eye(2);
I0 = zeros(2);

% Space
l1 = 0;             % x空間的下限
l2 = 1;             % x空間的上限
x_0_b = 0.25;
x_1_b = 0.75;
Delta{1} = 0.5;
Delta{2} = 0.5;
x_sample = 1/20;    % 將空間切成30份
zzz = (l1: x_sample: l2);  % Make the space interval. Different from before
N = round((l2-l1) / x_sample) + 1;

% Time
h_bar = 0.01;       % Maximum of h(t)

%% SOS program & Variable declear
g1 = 0.5;
g2 = 0.3;
rho{1} = 0.01;
rho{2} = 0.02;

po = sosprogram([y; r1; rr; rrr]);

[po, X]  = sospolymatrixvar(po, monomials([y], [1]), [2 2], 'symmetric');
for i = 1: d
    for j = 1: d
        if (i == j)
            X(i, j) = X(i, j);
        else
            X(i, j) = 0;
        end
    end 
end

[po, KB{1}{1}] = sospolymatrixvar(po, monomials([y], [2]), [1 2]);  % First number is v, second is j
[po, KB{1}{2}] = sospolymatrixvar(po, monomials([y], [2]), [1 2]);
[po, KB{2}{1}] = sospolymatrixvar(po, monomials([y], [2]), [1 2]);
[po, KB{2}{2}] = sospolymatrixvar(po, monomials([y], [2]), [1 2]);

[po, W1B] = sospolymatrixvar(po, monomials([y], [2]), [2 2], 'symmetric');
[po, W2B] = sospolymatrixvar(po, monomials([y], [2]), [2 2], 'symmetric');
[po, Q2B] = sospolymatrixvar(po, monomials([y], [2]), [2 2], 'symmetric');

[po, G] = sospolymatrixvar(po, monomials([y], [2]), [4 4], 'symmetric');
G11 = G(1:2, 1:2);
G12 = G(1:2, 3:4);
G21 = G(3:4, 1:2);
G22 = G(3:4, 3:4);

[po, Omega{1}] = sospolymatrixvar(po, monomials([y], [2]), [2 2], 'symmetric');
[po, Omega{2}] = sospolymatrixvar(po, monomials([y], [2]), [2 2], 'symmetric');


%% Condition
rn = 1;
% Positive defined matrix
w{rn} = rr' * X * rr;
po = sosineq(po, w{rn});
rn = rn + 1;
w{rn} = rr' * W1B * rr;
po = sosineq(po, w{rn});
rn = rn + 1;
w{rn} = rr' * W2B * rr;
po = sosineq(po, w{rn});
rn = rn + 1;
w{rn} = rr' * Q2B * rr;
po = sosineq(po, w{rn});
rn = rn + 1;
w{rn} = rr' * Omega{1} * rr;
po = sosineq(po, w{rn});
rn = rn + 1;
w{rn} = rr' * Omega{2} * rr;
po = sosineq(po, w{rn});
rn = rn + 1;

% % i = 1, j = 1, v = 0, v plus 1 可以從cell裡面取 ===========================
% i = 1;
% j = 1;
% v = 1;
% 
% PI11 = W1B + D(v)*KB{v}{j}*g2/Delta{v} + (D(v)*KB{v}{j}*g2/Delta{1})' + g2*A{i}*X + (g2*A{i}*X)' + rho{v}*Omega{v} - W2B;
% PI21 = Q2B + D{v}*KB{v}{j}*g1/Delta{v} - g2*X + g1*A{i}*X;
% PI22 = -2*g2*X;
% PI33 = -g2*Theta*X;
% PI41 = h_bar*G12;
% PI44 = -W1B + h_bar*G11;
% PI51 = h_bar*G12;
% PI55 = -h_bar*G22;
% PI61 = g2*KB{v}{j}'*D{v}'/Delta{v};
% PI62 = g1*KB{v}{j}'*D{v}'/Delta{v};
% PI66 = -Omega{v};
% PI71 = g2*KB{v}{j}'*D{v}'/Delta{v} + rho{v}*Omega{v};
% PI72 = g1*KB{v}{j}'*D{v}'/Delta{v};
% PI77 = -(pi^2) / (4*(Delta{v}^2)) * g2 * Theta * X + rho{v}*Omega{v};
% 
% PI = [PI11, PI21',   I0, PI41', PI51', PI61', PI71'; ...
%       PI21, PI22 ,   I0,   I0 ,   I0 , PI62', PI72'; ...
%         I0,   I0 , PI33,   I0 ,   I0 ,   I0 ,   I0 ; ...
%       PI41,   I0 ,   I0, PI44 ,   I0 ,   I0 ,   I0 ; ...
%       PI51,   I0 ,   I0,   I0 , PI55 ,   I0 ,   I0 ; ...
%       PI61, PI62 ,   I0,   I0 ,   I0 , PI66 ,   I0 ; ...
%       PI71, PI72 ,   I0,   I0 ,   I0 ,   I0 , PI77];
%   
% w{rn} = -r1' * PI * r1;
% po = sosineq(po, w{rn});
% rn = rn + 1;
% % =========================================================================
% 
% % i = 1, j = 1, v = 1, v plus 1 可以從cell裡面取 ===========================
% i = 1;
% j = 1;
% v = 2;
% 
% PI11 = W1B + D(v)*KB{v}{j}*g2/Delta{v} + (D(v)*KB{v}{j}*g2/Delta{1})' + g2*A{i}*X + (g2*A{i}*X)' + rho{v}*Omega{v} - W2B;
% PI21 = Q2B + D{v}*KB{v}{j}*g1/Delta{v} - g2*X + g1*A{i}*X;
% PI22 = -2*g2*X;
% PI33 = -g2*Theta*X;
% PI41 = h_bar*G12;
% PI44 = -W1B + h_bar*G11;
% PI51 = h_bar*G12;
% PI55 = -h_bar*G22;
% PI61 = g2*KB{v}{j}'*D{v}'/Delta{v};
% PI62 = g1*KB{v}{j}'*D{v}'/Delta{v};
% PI66 = -Omega{v};
% PI71 = g2*KB{v}{j}'*D{v}'/Delta{v} + rho{v}*Omega{v};
% PI72 = g1*KB{v}{j}'*D{v}'/Delta{v};
% PI77 = -(pi^2) / (4*(Delta{v}^2)) * g2 * Theta * X + rho{v}*Omega{v};
% 
% PI = [PI11, PI21',   I0, PI41', PI51', PI61', PI71'; ...
%       PI21, PI22 ,   I0,   I0 ,   I0 , PI62', PI72'; ...
%         I0,   I0 , PI33,   I0 ,   I0 ,   I0 ,   I0 ; ...
%       PI41,   I0 ,   I0, PI44 ,   I0 ,   I0 ,   I0 ; ...
%       PI51,   I0 ,   I0,   I0 , PI55 ,   I0 ,   I0 ; ...
%       PI61, PI62 ,   I0,   I0 ,   I0 , PI66 ,   I0 ; ...
%       PI71, PI72 ,   I0,   I0 ,   I0 ,   I0 , PI77];
%   
% w{rn} = -r1' * PI * r1;
% po = sosineq(po, w{rn});
% rn = rn + 1;
% % =========================================================================

% When hk = h_bar, h(t) = 0 ===============================================
for i = 1: 2
    for j = 1: 2
        for v = 1:2
            PI11 = W1B + D(v)*KB{v}{j}*g2/Delta{v} + (D(v)*KB{v}{j}*g2/Delta{1})' + g2*A{i}*X + (g2*A{i}*X)' + rho{v}*Omega{v} - W2B;
            PI21 = Q2B + D{v}*KB{v}{j}*g1/Delta{v} - g2*X + g1*A{i}*X;
            PI22 = -2*g2*X;
            PI33 = -g2*Theta*X;
            PI41 = h_bar*G12;
            PI44 = -W1B + h_bar*G11;
            PI51 = h_bar*G12;
            PI55 = -h_bar*G22;
            PI61 = g2*KB{v}{j}'*D{v}'/Delta{v};
            PI62 = g1*KB{v}{j}'*D{v}'/Delta{v};
            PI66 = -Omega{v};
            PI71 = g2*KB{v}{j}'*D{v}'/Delta{v} + rho{v}*Omega{v};
            PI72 = g1*KB{v}{j}'*D{v}'/Delta{v};
            PI77 = -(pi^2) / (4*(Delta{v}^2)) * g2 * Theta * X + rho{v}*Omega{v};

            PI = [PI11, PI21',   I0, PI41', PI51', PI61', PI71'; ...
                  PI21, PI22 ,   I0,   I0 ,   I0 , PI62', PI72'; ...
                    I0,   I0 , PI33,   I0 ,   I0 ,   I0 ,   I0 ; ...
                  PI41,   I0 ,   I0, PI44 ,   I0 ,   I0 ,   I0 ; ...
                  PI51,   I0 ,   I0,   I0 , PI55 ,   I0 ,   I0 ; ...
                  PI61, PI62 ,   I0,   I0 ,   I0 , PI66 ,   I0 ; ...
                  PI71, PI72 ,   I0,   I0 ,   I0 ,   I0 , PI77];
  
            w{rn} = -r1' * PI * r1;
            po = sosineq(po, w{rn});
            rn = rn + 1;     
        end
    end
end
% =========================================================================

% When hk = h_bar, h(t) = h_bar ===========================================
for i = 1: 2
    for j = 1: 2
        for v = 1:2
            PI11 = W1B + D(v)*KB{v}{j}*g2/Delta{v} + (D(v)*KB{v}{j}*g2/Delta{1})' + g2*A{i}*X + (g2*A{i}*X)' + rho{v}*Omega{v} - W2B;
            PI21 = Q2B + D{v}*KB{v}{j}*g1/Delta{v} - g2*X + g1*A{i}*X;
            PI22 = -2*g2*X;
            PI33 = -g2*Theta*X;
            PI41 = I0;
            PI44 = -W1B;
            PI51 = I0;
            PI55 = -h_bar*G22;
            PI61 = g2*KB{v}{j}'*D{v}'/Delta{v};
            PI62 = g1*KB{v}{j}'*D{v}'/Delta{v};
            PI66 = -Omega{v};
            PI71 = g2*KB{v}{j}'*D{v}'/Delta{v} + rho{v}*Omega{v};
            PI72 = g1*KB{v}{j}'*D{v}'/Delta{v};
            PI77 = -(pi^2) / (4*(Delta{v}^2)) * g2 * Theta * X + rho{v}*Omega{v};

            PI = [PI11, PI21',   I0, PI41', PI51', PI61', PI71'; ...
                  PI21, PI22 ,   I0,   I0 ,   I0 , PI62', PI72'; ...
                    I0,   I0 , PI33,   I0 ,   I0 ,   I0 ,   I0 ; ...
                  PI41,   I0 ,   I0, PI44 ,   I0 ,   I0 ,   I0 ; ...
                  PI51,   I0 ,   I0,   I0 , PI55 ,   I0 ,   I0 ; ...
                  PI61, PI62 ,   I0,   I0 ,   I0 , PI66 ,   I0 ; ...
                  PI71, PI72 ,   I0,   I0 ,   I0 ,   I0 , PI77];
  
            w{rn} = -r1' * PI * r1;
            po = sosineq(po, w{rn});
            rn = rn + 1;     
        end
    end
end
% =========================================================================

% When hk = 0, h(t) = 0. It leads to a reduced dimension. =================
for i = 1: 2
    for j = 1: 2
        for v = 1:2
            PI11 = W1B + D(v)*KB{v}{j}*g2/Delta{v} + (D(v)*KB{v}{j}*g2/Delta{1})' + g2*A{i}*X + (g2*A{i}*X)' + rho{v}*Omega{v} - W2B;
            PI21 = Q2B + D{v}*KB{v}{j}*g1/Delta{v} - g2*X + g1*A{i}*X;
            PI22 = -2*g2*X;
            PI33 = -g2*Theta*X;
            PI41 = I0;
            PI44 = -W1B;
            PI51 = I0;
            PI55 = I0;
            PI61 = g2*KB{v}{j}'*D{v}'/Delta{v};
            PI62 = g1*KB{v}{j}'*D{v}'/Delta{v};
            PI66 = -Omega{v};
            PI71 = g2*KB{v}{j}'*D{v}'/Delta{v} + rho{v}*Omega{v};
            PI72 = g1*KB{v}{j}'*D{v}'/Delta{v};
            PI77 = -(pi^2) * g2 / (4*(Delta{v}^2)) * Theta * X + rho{v}*Omega{v};

            PII = [PI11, PI21',  I0, PI41', PI61', PI71'; ...
                  PI21, PI22 ,   I0,   I0 , PI62', PI72'; ...
                    I0,   I0 , PI33,   I0 ,   I0 ,   I0 ; ...
                  PI41,   I0 ,   I0, PI44 ,   I0 ,   I0 ; ... ...
                  PI61, PI62 ,   I0,   I0 , PI66 ,   I0 ; ...
                  PI71, PI72 ,   I0,   I0 ,   I0 , PI77];
  
            w{rn} = -rrr' * PII * rrr;
            po = sosineq(po, w{rn});
            rn = rn + 1;     
        end
    end
end
% =========================================================================

save FHN_condition_finish.mat

%% Call solver
[po, info] = sossolve(po);

for j = 1: 2
    for v = 1: 2
        KB{v}{j}  = sosgetsol(po, KB{v}{j});
        KB{v}{j}  = sosgetsol(po, KB{v}{j});
    end
end

X = sosgetsol(po, X);
W1B = sosgetsol(po, W1B);
W2B = sosgetsol(po, W2B);
Q2B = sosgetsol(po, Q2B);
Omega{1} = sosgetsol(po, Omega{1});
Omega{2} = sosgetsol(po, Omega{2});


for j = 1:2
    for v = 1: 2
        K{v}{j} = KB{v}{j} * inv(X);
        K{v}{j} = KB{v}{j} * inv(X);
    end
end

save FHN_finish_calling_solver

%% Function
function [a] = a(g) % g is the position   
    d = 2;
    a = [zeros(d, d * (g-1)), eye(d), zeros(d, d * (7-g))]';
end
