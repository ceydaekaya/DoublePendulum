function doublePendulum(L1,L2,m1,m2,theta1_0,theta2_0)
L1=20;
L2=30;
theta1_0=pi/2;
theta2_0=3*pi/4;
m1=2;
m2=3;
tracing = 10;
g = 9.81;

% initial pose 
x = [mod(theta1_0,2*pi),mod(theta2_0,2*pi),0,0];
t = 0;
i_end = 0.02; %iteration duration
tracing = repmat(L1*[sin(x(1)),-cos(x(1))],[tracing,2]) + repmat([0 0 L2*[sin(x(2)),-cos(x(2))]],[tracing,1]);% trace of initial pose 

double_Pendulum = @(t,x)double_pendulum_system(x,L1,L2,m1,m2,g);% ODE function, i found this online 

axis xy equal, box on, hold on  
axis(1.1*[-1 1 -1 1]*(L1+L2))%to keep the pendulum inside the screen
r1=2;r2=2;%radius
iter = 0;

tic %starting the timer 
while 1 
    %solves differential equations for pendulum
    [t,x] = ode45(double_Pendulum,t(end)*[1 0.5 0] + i_end*[0 0.5 1] ,x(end,:)');
    tracing = patch_double_pendulum(t,x(end,:),L1,L2,r1,r2,tracing);%updating tracing to draw new points 
    iter = iter+1; %u for moving to next iteration
    i_end = max(toc*(1+1/iter),t(end)+2*eps); % end time for next iteration toc : checking how much time has elapsed
end
end

function dx = double_pendulum_system(x,L1,L2,m1,m2,g)% http://www.myphysicslab.com/dbl_pendulum.html
theta1 = x(1);
theta2 = x(2);
omega1 = x(3);%angular velocity of top rod
omega2 = x(4);%angular velocity of bottom rod
dtheta1 = omega1; %derivative of theta1 with respect to time
dtheta2 = omega2;
%equations of motion
domega1 = (-g*(2*m1+m2)*sin(theta1)-m2*g*sin(theta1-2*theta2)-2*sin(theta1-theta2)*m2*(omega2^2*L2+omega1^2*L1*cos(theta1-theta2)))/(L1*(2*m1+m2-m2*cos(2*theta1-2*theta2)));
domega2 = (2*sin(theta1-theta2)*(omega1^2*L1*(m1+m2)+g*(m1+m2)*cos(theta1)+omega2^2*L2*m2*cos(theta1-theta2)))/(L2*(2*m1+m2-m2*cos(2*theta1-2*theta2)));
dx = [dtheta1;dtheta2;domega1;domega2]; % diffential of x
end

function tracing = patch_double_pendulum(t,x,L1,L2,r1,r2,tracing)
cla %to prevent rode coordinates being traced
% Tracing
patch([tracing(:,1);NaN],[tracing(:,2);NaN],0,'EdgeColor','r','FaceColor','none','LineWidth',2);
patch([tracing(:,3);NaN],[tracing(:,4);NaN],0,'EdgeColor','g','FaceColor','none','LineWidth',2);
% plotting rods
theta1 = x(1);
theta2 = x(2);
xm1 = L1*sin(theta1);
ym1 = -L1*cos(theta1);
xm2 = xm1 + L2*sin(theta2);
ym2 = ym1 - L2*cos(theta2);
patch([0, xm1, xm2, NaN],[0, ym1, ym2, NaN],0,'EdgeColor','black','FaceColor','none','LineWidth',2)
% plotting bobs
p = linspace(0,2*pi,17);
sint = sin(p);
cost = cos(p);
patch(xm1+r1*cost,ym1+r1*sint,0,'EdgeColor','r','FaceColor','r')
patch(xm2+r2*cost,ym2+r2*sint,0,'EdgeColor','g','FaceColor','g')
drawnow
tracing = [tracing(1:end,:);xm1,ym1,xm2,ym2]; 
end

