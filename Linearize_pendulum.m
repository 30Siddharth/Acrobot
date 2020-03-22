init_point = [pi,0,0,0];

tspan = 0:0.001:1;
x0 = [3.145,0,0,0];
% x0 = init_point;
%% Linearizing the Model
syms q1 q2 q1dot q2dot 

% Parameters
m1 = 8;
m2 = 8;
l1 = 0.5;
l2 = 1;
g  = 9.81;
I1 = m1*l1^2;
I2 = m2*l2^2;


%%
M = [I1+I2+m2*l1^2 + 2*m2*l1*l2*cos(q2),  I2+ m2*l1*l2*cos(q2);
           I2 + m2*l1*l2*cos(q2),                      I2];
       
C = [-2*m2*l1*l2*sin(q2)*q2dot,  -m2*l1*l2*sin(q2)*q2dot;
     m2*l1*l2*sin(q2)*q1dot,              0];
 
T = [-m1*g*l1*sin(q1)-m2*g*(l1*sin(q1) + l2*sin(q1+q2));
            -m2*g*l2*sin(q1+q2)];
       
f2 = inv(M)*{T - C*[q1dot;q2dot]};

f1 = [q1dot;q2dot];

f = [f1;f2];
f0 = double(subs(f,[q1,q2,q1dot,q2dot], init_point));

%%
diff_f = jacobian(f,[q1,q2,q1dot,q2dot]);
A = double(subs(diff_f,[q1,q2,q1dot,q2dot], init_point));


%%
del_f_del_u = [0;
     0;
     inv(M)*[0;1]];
B = double(subs(del_f_del_u, [q1,q2,q1dot,q2dot], init_point));


% x


%% ODE
[t,state] = ode45(@(t,x) dynamics(t,x,A,B,f0,init_point),tspan,x0);
x1 = l1*sin(state(:,1));
y1 = -l1*cos(state(:,1));

x2 = x1 + l2*sin(state(:,1) + state(:,2));
y2 = y1 - l2*cos(state(:,1) + state(:,2));

figure(1)
plot(t,state,'LineWidth',3);
legend('q1','q2','q1dot','q2dot');
grid on


% figure(2)
% hold on;
% plot(x1,y1,'-*','MarkerIndices',1:100:length(y1));
% 
% % Plot of pendulum in starting postion
% plot([0,x1(1)],[0,y1(1)], 'LineWidth', 4)
% plot([x1(1),x2(1)],[y1(1),y2(1)], 'LineWidth', 4)
% %Plot of pendulum in end postion
% plot([0,x1(end)],[0,y1(end)], 'LineWidth', 4)
% plot([x1(end),x2(end)],[y1(end),y2(end)], 'LineWidth', 4)


% plot(x1(1),y1(1),'o','MarkerFaceColor','red','MarkerSize',6);
% plot(x1(end),y1(end),'o','MarkerFaceColor','red','MarkerSize',6);
% % 'MarkerFaceColor','red',...
%     
% figure(2)
% hold on;
% plot(x2,y2);
% plot(x2(1),y2(1),'o','MarkerFaceColor','red','MarkerSize',6);
% plot(x2(end),y2(end),'o','MarkerFaceColor','red','MarkerSize',6);


% plot(0,0,'-p','MarkerFaceColor','green','MarkerSize',10);
% grid on
% axis equal

%% Animation
k = length(tspan);
% for t = 1:50:k
%     
%     % Plot of pendulum in starting postion
%     figure(2)
%     hold on;
%     plot([0,x1(t)],[0,y1(t)], 'LineWidth', 4)
%     plot([x1(t),x2(t)],[y1(t),y2(t)], 'LineWidth', 4)
%     plot(x1(t),y1(t),'o','MarkerFaceColor','red','MarkerSize',6);
%     plot(x2(t),y2(t),'o','MarkerFaceColor','red','MarkerSize',6);
%     grid on
%     axis equal
%     
% end
%% Dynamics 
function dxdt = dynamics(t,x,A,B,f0,x0)
Q = 100*eye(4);
R = 1;

[X,K,L] = icare(A,B,Q,R);
u = -K*(x-x0');
x-x0'
dxdt = A*(x-x0') + B*u +f0;

end

        
