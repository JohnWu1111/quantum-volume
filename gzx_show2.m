clear;
clc;
format long
tic;

%% paramter

h0 = 0:0.01:2;
h1 = 0:0.01:2;
nh = length(h0);
dk = 2*pi/1000;
k = 0:dk:pi;
nk = length(k);

Gkk = 2*(sin(k')*h).^2./(1+h.^2+2*cos(k')*h);
QV = 2*sum(Gkk).*dk;

figure
mesh(h,k,Gkk)
xlabel('h')
ylabel('k')
% plot(h,QV)

toc;