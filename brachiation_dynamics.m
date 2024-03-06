% Script for symbolically compuing the dynamics of a 5DOF brachiating robot

m=1; % Mass of each link. All links are assumed to have the same mass
l=1; % Length of each link. All links are assumed to have the same length
g=9.8; % Gravity in m/s^2

% Configuration variables of the robot. Angles are measured
% counter-clockwise between adjacent links.
syms q1 q2 q3 q4 q5 q1dot q2dot q3dot q4dot q5dot real

q=[q1;q2;q3;q4;q5];
qdot=[q1dot;q2dot;q3dot;q4dot;q5dot];

% Centres of mass of each link. Each link is assumed to be a point mass of
% mass m attached to the end of a massless rod of length l
rc1=l*[-sin(q1);cos(q1)];
rc2=rc1+l*[-sin(q1+q2);cos(q1+q2)];
rc3=rc2+l*[-sin(q1+q2+q3);cos(q1+q2+q3)];
rc4=rc3+l*[-sin(q1+q2+q3+q4);cos(q1+q2+q3+q4)];
rc5=rc2+l*[-sin(q1+q2+q5);cos(q1+q2+q5)];

% Time derivatives of centres of mass
rc1dot=jacobian(rc1,q)*qdot;
rc2dot=jacobian(rc2,q)*qdot;
rc3dot=jacobian(rc3,q)*qdot;
rc4dot=jacobian(rc4,q)*qdot;
rc5dot=jacobian(rc5,q)*qdot;

rc=[rc1;rc2;rc3;rc4;rc5];
rcdot=[rc1dot;rc2dot;rc3dot;rc4dot;rc5dot];

% Total kinetic energy of the robot
K=1/2*m*(rc1dot'*rc1dot ...
           +rc2dot'*rc2dot ...
           +rc3dot'*rc3dot ...
           +rc4dot'*rc4dot ...
           +rc5dot'*rc5dot);
K=simplify(K)

D=simplify(hessian(K,qdot))
C = sym(zeros(5,5));
for i=1:5
    for j=1:5
        for k=1:5
            C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qdot(i);
        end
    end
end

% Total potential energy of the robot
P = m*g*(rc1(2)+rc2(2)+rc3(2)+rc4(2)+rc5(2))

dP = jacobian(P,q)'

% Input matrix
B = [0,0,0,0;
     1,0,0,0;
     0,1,0,0;
     0,0,1,0;
     0,0,0,1];

% Troques appplied at each joint of the robot
syms u1 u2 u3 u4 real
u = [u1;u2;u3;u4];

%System dynamics
ddq = D\(B*u-C*qdot-dP)
