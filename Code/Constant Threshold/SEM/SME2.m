clc; clear; close all;
tic; echo off;
% load SEM_finish_calling_solver.mat
load test7.mat
% load ada0503_4.mat

%% Time & Space 
t(1) = 0;
t0   = 0;    
tf   = 40;
t_sample = 0.00125;
N_t = round((tf-t0) / t_sample) + 1;  % total step
ttt = t0: t_sample: tf;

l1 = 0;
l2 = 1;
xp1 = 5;
xp2 = 16;
zzz(1) = 0;

Y1(1,1)  = 0;
Y2(1,1)  = 0;
yy1(1,1) = 0;
yy2(1,1) = 0;

Y1  = repmat(Y1, N_t, 1);
Y2  = repmat(Y2, N_t, 1);
yy1 = repmat(yy1, N, 1);
yy2 = repmat(yy2, N, 1);
u1  = zeros(N_t, 2);
u2  = zeros(N_t, 2);
rhos1 = zeros(N_t, 1);
rhos2 = zeros(N_t, 1);
rho1 = zeros(N_t, 1);
rho1(1) = 0.06;
rho2 = zeros(N_t, 1);
rho2(1) = 0.07;
rc = 0.01;

%% Initial condition
for i = 1: N
    Y1(1,i)  =  0.5 * sin(pi*zzz(i));
    Y2(1,i)  =  0.3 * sin(pi*zzz(i));        %  0.1 * sin(pi*zzz(i)) at first
    yy1(1,i) = -0.5 * pi^2 * sin(pi*zzz(i));
    yy2(1,i) = -0.3 * pi^2 * sin(pi*zzz(i)); %  0.1 * sin(pi*zzz(i)) at first
end

% The selected states for controller 
y1tkv1 = Y1(1,xp1); 
y2tkv1 = Y2(1,xp1);
y1tkv2 = Y1(1,xp2);
y2tkv2 = Y2(1,xp2);

%% Dynamic 
for it = 1: N_t-1
    for in = 1: N
        % Mechanism at v = 0
%         if( rho1(it) >= rho{1} )
%             % Over the upper bound
%             rho1(it) = rho{1};
%         else
%             if ( norm([y1tkv1; y2tkv1]) <= norm([Y1(it,xp1); Y2(it,xp1)]) )
%                 % Lower bound
%                 rho1(it) = rho1(1);
%             else
%                 % State converge => threshold increasing
%                 nl = tanh( (norm([Y1(it,xp1); Y2(it,xp1)]) - norm([y1tkv1; y2tkv1])) /  norm([Y1(it,xp1); Y2(it,xp1)]) );
%                 rho1(it) = (1 - rc*nl)*rho1(it-1);
%             end
%         end               
        
        if ( (it>=2) && (in == xp1) && ([y1tkv1-Y1(it, xp1); y2tkv1-Y2(it, xp1)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp1); Y2(it, xp1)]) *  [y1tkv1-Y1(it, xp1); y2tkv1-Y2(it, xp1)] < rho{1} * [Y1(it, xp1); Y2(it, xp1)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp1); Y2(it, xp1)]) *  [Y1(it, xp1); Y2(it, xp1)]))
            fprintf("Threshold : %d. Does not pull the trigger at v1 (%d. %d). \n", rho1(it), it, in)
        elseif ( (it>=2) && (in == xp1) && ([y1tkv1-Y1(it, xp1); y2tkv1-Y2(it, xp1)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp1); Y2(it, xp1)]) *  [y1tkv1-Y1(it, xp1); y2tkv1-Y2(it, xp1)] >= rho{1} * [Y1(it, xp1); Y2(it, xp1)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp1); Y2(it, xp1)]) *  [Y1(it, xp1); Y2(it, xp1)]))
            fprintf("Threshold : %d. Pull the trigger at v1 (%d. %d). \n", rho1(it), it, in)
            y1tkv1 = Y1(it, xp1);
            y2tkv1 = Y2(it, xp1);
            rhos1(it) = 1;
        end
        % Mechanism at v = 1
%         if( rho2(it) >= rho{2} )
%             % Over the upper bound
%             rho2(it) = rho{2};
%         else
%             if ( norm([y1tkv2; y2tkv2]) <= norm([Y1(it,xp2); Y2(it,xp2)]) )
%                 % Lower bound
%                 rho2(it) = rho2(1);
%             else
%                 % State converge => threshold increasing
%                 nl = tanh( (norm([Y1(it,xp2); Y2(it,xp2)]) - norm([y1tkv2; y2tkv2])) /  norm([Y1(it,xp2); Y2(it,xp2)]) );
%                 rho2(it) = (1 - rc*nl)*rho2(it-1);
%             end
%         end
        
        if ( (it>=2) && (in == xp2) && ([y1tkv2-Y1(it, xp2); y2tkv2-Y2(it, xp2)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp2); Y2(it, xp2)]) *  [y1tkv2-Y1(it, xp2); y2tkv2-Y2(it, xp2)] < rho{2} *[Y1(it, xp2); Y2(it, xp2)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp2); Y2(it, xp2)]) *  [Y1(it, xp2); Y2(it, xp2)]))
            fprintf("Threshold : %d. Does not pull the trigger at v2 (%d. %d). \n", rho2(it), it, in)
        elseif ( (it>=2) && (in == xp2) && ([y1tkv2-Y1(it, xp2); y2tkv2-Y2(it, xp2)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp2); Y2(it, xp2)]) *  [y1tkv2-Y1(it, xp2); y2tkv2-Y2(it, xp2)] >= rho{2} *[Y1(it, xp2); Y2(it, xp2)]' * subs(Omega{v}, [y1; y2], [Y1(it, xp2); Y2(it, xp2)]) *  [Y1(it, xp2); Y2(it, xp2)]))
            fprintf("Threshold : %d. Pull the trigger at v2 (%d. %d). \n", rho2(it), it, in)
            y1tkv2 = Y1(it, xp2);
            y2tkv2 = Y2(it, xp2);
            rhos2(it) = 1;
        end
        
        % Membership function
        h1 = (sin(Y1(it, in)) - Y1(it, in) + 0.1667*Y1(it, in)^3) / 0.0178*Y1(it, in)^3;
        h2 = 1 - h1;
        % Gain the double form
        KF11 = double(subs(K{1}{1}, [y1, y2], [y1tkv1, y2tkv1]));
        KF12 = double(subs(K{1}{2}, [y1, y2], [y1tkv1, y2tkv1]));
        KF21 = double(subs(K{2}{1}, [y1, y2], [y1tkv2, y2tkv2]));
        KF22 = double(subs(K{2}{2}, [y1, y2], [y1tkv2, y2tkv2]));    
        
        u{1} = h1*KF11*[y1tkv1; y2tkv1] + h2*KF12*[y1tkv1; y2tkv1]; % u at v1
        u{2} = h1*KF21*[y1tkv2; y2tkv2] + h2*KF22*[y1tkv2; y2tkv2];
        u1(it, :) = u{1};
        u2(it, :) = u{2};
        
        % Euler
        if (in*x_sample <= 0.5)
            v = 1;
            Y1(it+1, in) = Y1(it, in) + t_sample * subs( Theta(1, 1)*yy1(it,in) + h1*A{1}(1,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(1,:)*[Y1(it,in); Y2(it,in)] + D{v}(1,:)*u{v}, [y1; y2], [Y1(it,in); Y2(it,in)]  );
            Y2(it+1, in) = Y2(it, in) + t_sample * subs( Theta(2, 2)*yy2(it,in) + h1*A{1}(2,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(2,:)*[Y1(it,in); Y2(it,in)] + D{v}(2,:)*u{v}, [y1; y2], [Y1(it,in); Y2(it,in)]  );
            
        else 
            v = 2;
            Y1(it+1, in) = Y1(it, in) + t_sample * subs( Theta(1, 1)*yy1(it,in) + h1*A{1}(1,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(1,:)*[Y1(it,in); Y2(it,in)] + D{v}(1,:)*u{v}, [y1; y2], [Y1(it,in); Y2(it,in)]  );
            Y2(it+1, in) = Y2(it, in) + t_sample * subs( Theta(2, 2)*yy2(it,in) + h1*A{1}(2,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(2,:)*[Y1(it,in); Y2(it,in)] + D{v}(2,:)*u{v}, [y1; y2], [Y1(it,in); Y2(it,in)]  );
            
        end
        
    end % for in
    
    for i = 2: N-1
        yy1(it+1,i) = (Y1(it+1,i+1) - 2*Y1(it+1,i) + Y1(it+1,i-1)) / x_sample / x_sample;
        yy2(it+1,i) = (Y2(it+1,i+1) - 2*Y2(it+1,i) + Y2(it+1,i-1)) / x_sample / x_sample;           
    end
    
    
end % for it

%% Figure
figure
set(gcf, 'Renderer', 'ZBuffer');
mesh(zzz, ttt, Y1)
view(-40+90, 30);
xlabel('x');
ylabel('t');
zlabel('y_1');

figure
set(gcf, 'Renderer', 'ZBuffer');
mesh(zzz, ttt, Y2)
view(-40+90, 30);
xlabel('x');
ylabel('t');
zlabel('y_2');

figure
plot(ttt, rho1); hold on
plot(ttt, rho2);
legend("rho v1", "rho v2")

figure
plot(ttt, u1(:,1), ttt, u1(:,2), ttt, u2(:,1), ttt, u2(:,2)); hold on;
legend("u11", "u12", "u21", "u22");




toc
