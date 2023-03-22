clc; clear; close all;

load FHN_finish_calling_solver.mat

%% Time & Space 
t(1) = 0;
t0   = 0;    
tf   = 20;
t_sample = 0.001;
N_t = round((tf-t0) / t_sample) + 1;  % total step
ttt = t0: t_sample: tf;

l1 = 0;
l2 = 1;
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

%% Dynamic 
for it = 1: N_t-1
    for in = 1: N
        h1 = alpha^(-2) * Y1(it, in);
        h2 = 1 - h1;
        
        if (in*x_sample <= 0.5)
            v = 1;
            Y1(it+1, in) = Y1(it, in) + subs( t_sample * ( Theta(1, 1)*yy1(it,in) + h1*h1*(A{1}(1,:)+D{v}(1)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h1*h2*(A{1}(1,:)+D{v}(1)*K{v}{2})*[Y1(it,in); Y2(it,in)] + h2*h1*(A{2}(1,:)+D{v}(1)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h2*h2*(A{2}(1,:)+D{v}(1)*K{v}{2})*[Y1(it,in); Y2(it,in)] ), [y1; y2], [Y1(it,6); Y2(it,6)]);
            Y2(it+1, in) = Y2(it, in) + subs( t_sample * ( Theta(2, 2)*yy2(it,in) + h1*h1*(A{1}(2,:)+D{v}(2)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h1*h2*(A{1}(2,:)+D{v}(2)*K{v}{2})*[Y1(it,in); Y2(it,in)] + h2*h1*(A{2}(2,:)+D{v}(2)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h2*h2*(A{2}(2,:)+D{v}(2)*K{v}{2})*[Y1(it,in); Y2(it,in)] ), [y1; y2], [Y1(it,6); Y2(it,6)]);
            
        else 
            v = 2;
            Y1(it+1, in) = Y1(it, in) + subs( t_sample * ( Theta(1, 1)*yy1(it,in) + h1*h1*(A{1}(1,:)+D{v}(1)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h1*h2*(A{1}(1,:)+D{v}(1)*K{v}{2})*[Y1(it,in); Y2(it,in)] + h2*h1*(A{2}(1,:)+D{v}(1)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h2*h2*(A{2}(1,:)+D{v}(1)*K{v}{2})*[Y1(it,in); Y2(it,in)] ), [y1; y2], [Y1(it,6); Y2(it,16)]);
            Y2(it+1, in) = Y2(it, in) + subs( t_sample * ( Theta(2, 2)*yy2(it,in) + h1*h1*(A{1}(2,:)+D{v}(2)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h1*h2*(A{1}(2,:)+D{v}(2)*K{v}{2})*[Y1(it,in); Y2(it,in)] + h2*h1*(A{2}(2,:)+D{v}(2)*K{v}{1})*[Y1(it,in); Y2(it,in)] + h2*h2*(A{2}(2,:)+D{v}(2)*K{v}{2})*[Y1(it,in); Y2(it,in)] ), [y1; y2], [Y1(it,6); Y2(it,16)]);
            
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



