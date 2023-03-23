clc; clear; close all;

load FHN_finish_calling_solver.mat

%% Time & Space 
t(1) = 0;
t0   = 0;    
tf   = 20;
t_sample = 0.00125;
N_t = round((tf-t0) / t_sample) + 1;  % total step
ttt = t0: t_sample: tf;

l1 = 0;
l2 = 1;
xp1 = 6;
xp2 = 16;
x_sample = 1 / 20;
zzz(1) = 0;

Y1(1,1)  = 0;
Y2(1,1)  = 0;
yy1(1,1) = 0;
yy2(1,1) = 0;

Y1  = repmat(Y1, N_t, 1);
Y2  = repmat(Y2, N_t, 1);
yy1 = repmat(yy1,  N, 1);
yy2 = repmat(yy2,  N, 1);
u1  = zeros(N_t,1);
u2  = zeros(N_t,1);

%% Initial condition
for i = 1: N
    Y1(1,i)  =  0.5 * cos(pi*zzz(i)) + 0.5;
    Y2(1,i)  =  0.1 * cos(pi*zzz(i));
    yy1(1,i) = -0.5 * pi^2 * cos(pi*zzz(i));
    yy2(1,i) = -0.1 * pi^2 * cos(pi*zzz(i));
end

% The selected states for controller 
y1tkv1 = Y1(1,xp1); 
y2tkv1 = Y2(1,xp1);
y1tkv2 = Y1(1,xp2);
y2tkv2 = Y2(1,xp2);

%% Dynamic 
for it = 1: N_t-1
    for in = 1: N
        h1 = alpha^(-2) * Y1(it, in);
        h2 = 1 - h1;
        
        % Mechanism at v = 0
        if ( (it>=2) && (in == xp1) && [y1tkv1-Y1(it, xp1); y2tkv1-Y2(it, xp1)]' * Omega{1} *  [y1tkv1-Y1(it, xp1); y2tkv1-Y2(it, xp1)] < rho{1} *[Y1(it, xp1); Y2(it, xp1)]' * Omega{1} *  [Y1(it, xp1); Y2(it, xp1)])
            fprintf("Does not pull the triggered at v1")
        else
            y1tkv1 = Y1(it, xp1);
            y2tkv1 = Y2(it, xp1);
        end
        % Mechanism at v = 1
        if ( (it>=2) && (in == xp2) && [y1tkv1-Y1(it, xp2); y2thv1-Y2(it, xp2)]' * Omega{1} *  [y1tkv1-Y1(it, xp2); y2thv1-Y2(it, xp2)] < rho{1} *[Y1(it, xp2); Y2(it, xp2)]' * Omega{1} *  [Y1(it, xp2); Y2(it, xp2)])
            fprintf("Does not pull the triggered at v2")
        else
            y1tkv2 = Y1(it, xp2);
            y2tkv2 = Y2(it, xp2);
        end
        
        % Euler
        if (in*x_sample <= 0.5)
            v = 1;
            Y1(it+1, in) = Y1(it, in) + subs( t_sample * ( Theta(1, 1)*yy1(it,in) + h1*A{1}(1,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(1,:)*[Y1(it,in); Y2(it,in)] + D{v}(1,:)*(h1*K{v}{1}*[y1tkv1; y2tkv1] +  h2*K{v}{2}*[y1tkv1; y2tkv1]) ), [y1; y2], [y1tkv1; y2tkv1]);
            Y2(it+1, in) = Y2(it, in) + subs( t_sample * ( Theta(2, 2)*yy2(it,in) + h1*A{1}(2,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(2,:)*[Y1(it,in); Y2(it,in)] + D{v}(2,:)*(h1*K{v}{1}*[y1tkv1; y2tkv1] +  h2*K{v}{2}*[y1tkv1; y2tkv1]) ), [y1; y2], [y1tkv1; y2tkv1]);
            
        else 
            v = 2;
            Y1(it+1, in) = Y1(it, in) + subs( t_sample * ( Theta(1, 1)*yy1(it,in) + h1*A{1}(1,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(1,:)*[Y1(it,in); Y2(it,in)] + D{v}(1,:)*(h1*K{v}{1}*[y1tkv1; y2tkv1] +  h2*K{v}{2}*[y1tkv1; y2tkv1]) ), [y1; y2], [y1tkv1; y2tkv1]);
            Y2(it+1, in) = Y2(it, in) + subs( t_sample * ( Theta(2, 2)*yy2(it,in) + h1*A{1}(2,:)*[Y1(it,in); Y2(it,in)] + h2*A{2}(2,:)*[Y1(it,in); Y2(it,in)] + D{v}(2,:)*(h1*K{v}{1}*[y1tkv1; y2tkv1] +  h2*K{v}{2}*[y1tkv1; y2tkv1]) ), [y1; y2], [y1tkv1; y2tkv1]);
            
        end
        
    end % for in
    
    for i = 2: N-1
        if i == 1
            yy1(it+1,i) = 0;
            yy2(it+1,i) = 0;
        else
            yy1(it+1,i) = (Y1(it+1,i+1) - 2*Y1(it+1,i) + Y1(it+1,i-1)) / x_sample / x_sample;
            yy2(it+1,i) = (Y2(it+1,i+1) - 2*Y2(it+1,i) + Y2(it+1,i-1)) / x_sample / x_sample;           
        end
    end
    
    
end % for it

%% Figure
figure(1)
set(gcf, 'Renderer', 'ZBuffer');
mesh(zzz, ttt, Y1)
view(-40+90, 30);
xlabel('x');
ylabel('t');
zlabel('y_1');

figure(2)
set(gcf, 'Renderer', 'ZBuffer');
mesh(zzz, ttt, Y2)
view(-40+90, 30);
xlabel('x');
ylabel('t');
zlabel('y_2');



