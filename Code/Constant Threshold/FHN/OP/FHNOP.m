clear; clc; close all;

tic

l1 = 0;
l2 = 1;
x_sample = 1 / 15;
zzz(1) = 0;

alpha = 1; % alpha = 0.0001
N = round((l2-l1) / x_sample) + 1;

for i = 1: N
    zzz(i) = (i-1) * x_sample;
end

t0 = 0;
tf = 30;
t_sample = 0.00125; % Normal: 0.001 sec. Max test: 0.00125 

N_t = round((tf-t0) / t_sample) + 1;
ttt = ((1:N_t)-1) * t_sample;

Y1(1,1) = 0;
Y2(1,1) = 0;
yy1(1,1) = 0;
yy2(1,1) = 0;
% yyy1(1,1) = 0;
% yyy2(1,1) = 0;
Y1 = repmat(Y1, N_t, 1);
Y2 = repmat(Y2, N_t, 1);
yy1 = repmat(yy1, N, 1);
yy2 = repmat(yy2, N, 1);
% yyy1 = repmat(yy1, N, 1);
% yyy2 = repmat(yy2, N, 1);



for i = 1: N
    % Initial condition
    Y1(1,i)  =  0.5 * cos(pi*zzz(i)) + 0.5;
    Y2(1,i)  = -0.3 * cos(pi*zzz(i));
    yy1(1,i) = -0.5 * pi^2 * cos(pi*zzz(i));
    yy2(1,i) =  0.3 * pi^2 * cos(pi*zzz(i));
end

for it = 1: (N_t-1)
    for i = 1: N
        Y1(it+1, i) = Y1(it,i) + t_sample * (alpha*yy1(it,i) + Y1(it,i) - Y2(it,i) - (Y1(it,i)^3));
        Y2(it+1, i) = Y2(it,i) + t_sample * (alpha*yy2(it,i) + 0.45*Y1(it,i) - 0.1*Y2(it,i));        
    end 
    
    for i = 2: N-1
        if i == 1
            yy1(it+1,i) = 0;
            yy2(it+1,i) = 0;
        else
            yy1(it+1,i) = (Y1(it+1,i+1) - 2*Y1(it+1,i) + Y1(it+1,i-1)) / x_sample / x_sample;
            yy2(it+1,i) = (Y2(it+1,i+1) - 2*Y2(it+1,i) + Y2(it+1,i-1)) / x_sample / x_sample;           
        end
    end   
end

figure
set(gcf, 'Renderer', 'ZBuffer');
mesh(zzz, ttt, Y1)
view(-40+90, 30);
xlabel('$x$', 'Interpreter','latex');
ylabel('$t$', 'Interpreter','latex');
zlabel('$y_1$', 'Interpreter','latex');

figure
set(gcf, 'Renderer', 'ZBuffer');
mesh(zzz, ttt, Y2)
view(-40+90, 30);
xlabel('$x$', 'Interpreter','latex');
ylabel('$t$', 'Interpreter','latex');
zlabel('$y_2$', 'Interpreter','latex');

save FHNOP.mat

toc


% tic
% 
% l1 = 0;
% l2 = 1;
% x_sample = 1 / 29;
% zzz(1) = 0;
% 
% alpha = -0.1;
% N = round((l2-l1) / x_sample) + 1;
% 
% for i = 1: N
%     zzz(i) = (i-1) * x_sample;
% end
% 
% t0 = 0;
% tf = 30;
% t_sample = 0.001; % First t:0.001 Second t:?
% 
% N_t = round((tf-t0) / t_sample) + 1;
% ttt = ((1:N_t)-1) * t_sample;
% 
% Y1(1,1) = 0;
% Y2(1,1) = 0;
% yy1(1,1) = 0;
% yy2(1,1) = 0;
% Y1 = repmat(Y1, N_t, 1);
% Y2 = repmat(Y2, N_t, 1);
% yy1 = repmat(yy1, N, 1);
% yy2 = repmat(yy2, N, 1);
% 
% 
% for i = 1: N
%      % Initial condition
%      Y1(1,i)   = -0.4 * sin(pi*zzz(i));
%      Y2(1,i)   =  0.2 * sin(pi*zzz(i));
%      yy1(1,i) =  0.4 * pi^2 * sin(pi*zzz(i));
%      yy2(1,i) = -0.2 * pi^2 * cos(pi*zzz(i));
% end
% 
%  for it = 1: (N_t-1)
%      for i = 1: N
%          Y1(it+1, i) = Y1(it,i) + t_sample * (alpha*yy1(it,i) + 3*Y1(it,i) - Y2(it,i) - (Y1(it,i)^3));
%          Y2(it+1, i) = Y2(it,i) + t_sample * (alpha*yy2(it,i) + 0.45*Y1(it,i) - 0.1*Y2(it,i));        
%      
%          if i == 1
%             yy1(it+1,i) = 0;
%             yy2(it+1,i) = 0;
%          else
%             yy1(it+1,i) = (Y1(it+1,i) - Y1(it+1,i-1)) / x_sample;
%             yy2(it+1,i) = (Y2(it+1,i) - Y2(it+1,i-1)) / x_sample;
%          end 
%      end
% end
% 
% figure(1)
% %set(gcf, 'Renderer', 'ZBuffer');
% mesh(zzz, ttt, Y1)
% view(-40+90, 30);
% xlabel('x');
% ylabel('t');
% zlabel('y_1');
% 
% figure(2)
% set(gcf, 'Renderer', 'ZBuffer');
% mesh(zzz, ttt, Y2)
% view(-40+90, 30);
% xlabel('x');
% ylabel('t');
% zlabel('y_2');
%  
% toc

