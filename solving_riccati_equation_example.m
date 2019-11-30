% This file serves as an example to program to solve riccati equation
% The implementation is based on the example explained by Christopher Lum
% from his video "Introduction to Linear Quadratic Regulator (LQR) Control"
% link: https://youtu.be/wEevt2a4SKI
clc
clear all
close all
%% Step 1: Define A and B
A=[0 1;0 -1/5] % note that A is a singular matrix here
B=[0;1]

%% Step 2: Choose Q and R
Q=[1 0;0 1]
R=[1/100]

%% Step 3: Solve Algebraic Riccati equation
% A'S + SA - SBinv(R)B'S + Q = 0
syms s11 s12 s22
% Define S as symmetric matrix
S=[s11 s12; s12 s22]
LHS=A'*S+S*A-S*B*inv(R)*B'*S+Q
% Since RHS = 0 then each of the element is 0 
sol=solve(LHS(1,1)==0,LHS(1,2)==0,LHS(2,2)==0,s11,s12,s22)
% we have 4 solutions to the algebraic riccati equation
S1=[sol.s11(1) sol.s12(1);sol.s12(1) sol.s22(1)];
S1=double(S1)
S2=[sol.s11(2) sol.s12(2);sol.s12(2) sol.s22(2)];
S2=double(S2)
S3=[sol.s11(3) sol.s12(3);sol.s12(3) sol.s22(3)];
S3=double(S3)
S4=[sol.s11(4) sol.s12(4);sol.s12(4) sol.s22(4)];
S4=double(S4)

%% Step 4: Compute K
% K = inv(R)B'S
K1 = inv(R)*B'*S1
K2 = inv(R)*B'*S2
K3 = inv(R)*B'*S3
K4 = inv(R)*B'*S4

%% Step 5: Find solution that yields stable system
E1=eig(A-B*K1)
E2=eig(A-B*K2)
E3=eig(A-B*K3)
E4=eig(A-B*K4) % stable solution
% Cross checking with in-built MATLAB command
[K_m,S_m,E_m]=lqr(A,B,Q,R)
% and it matches!