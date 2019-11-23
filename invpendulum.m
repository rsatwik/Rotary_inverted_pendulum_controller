clc
clear all
close all
Mp = 0.027
lp = 0.153
Lp = 0.191
r = 0.08260
Jm = 3e-5
Marm = 0.028
g = 9.810
Jeq = 1.23e-4
Jp = 1.1e-4
Beq = 0
Bp = 0
Rm = 3.30
Kt = 0.02797
Km = 0.02797
A = [0 0 1 0;0 0 0 1;0 (r*(Mp^2)*(lp^2)*g)/(Jp*Jeq + Mp*(lp^2)*Jeq + Jp*Mp*(r^2)) -(Kt*Km)*(Jp+Mp*(lp^2))/((Jp*Jeq + Mp*(lp^2)*Jeq + Jp*Mp*(r^2))*(Rm)) 0;0 (Mp*lp*g*(Jeq + Mp*(r^2)))/(Jp*Jeq + Mp*(lp^2)*Jeq + Jp*Mp*(r^2)) -1*(Mp*lp*Kt*r*Km)/((Jp*Jeq + Mp*(lp^2)*Jeq + Jp*Mp*(r^2))*Rm) 0]
B = [0;0;Kt*(Jp+Mp*(lp^2))/((Jp*Jeq + Mp*(lp^2)*Jeq + Jp*Mp*(r^2))*Rm);(Mp*lp*Kt*r)/((Jp*Jeq + Mp*(lp^2)*Jeq + Jp*Mp*(r^2))*Rm)]
C=eye(4)
D=zeros(4,1)
Q=diag([1,0,0,0])
R=1
format short;
K=lqr(A,B,Q,R)
C1=[1 0 0 0];
D1=[0];
%sys2=ss((A-B*K),[0;0;0;0],C1,D1);
%initial(sys2,[0;0.052;0;0])
%legend(['Q=' mat2str(Q)],['R=' num2str(R)]);
%%
%hold on
Q=diag([500 700 20 20]);
R=11;
K=lqr(A,B,Q,R)
%eig(A-B*K)
%C1=[1 0 0 0];
%D1=[0];
%sys2=ss((A-B*K),[0;0;0;0],C1,D1);
%initial(sys2,[0;0.052;0;0])
%legend(['Q=' mat2str(Q)],['R=' num2str(R)]);
%%
C1=[1 0 0 0];
D1=[0];
sys2=ss(A,B,C1,D1);
initial(sys2,[0;0.5;0;0])