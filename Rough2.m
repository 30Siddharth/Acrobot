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
       
% g = [0;
%     0;
%     inv(M)*[0;1]];
% x0 = [0,0,0,0];
% B = double(subs(g, [q1,q2,q1dot,q2dot],x0))
del_f_del_u = [0;
     0;
     inv(M)*[0;1]]
% B = double(subs(del_f_del_u, [q1,q2,q1dot,q2dot], init_point));