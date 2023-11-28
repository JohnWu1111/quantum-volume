clear;
clc;
format long
tic;

L = 1e5;
dk = 2*pi/L;
k = 0:dk:2*pi-dk;

phi = randn(L,2);
phi = phi./sqrt(abs(phi(:,1)).^2+abs(phi(:,2)).^2);

dphi = (circshift(phi,-1,1) - circshift(phi,1,1))/(2*dk);

tar = sum(sqrt(sum(dphi.*dphi - 2*dphi.*phi,2)));

toc;