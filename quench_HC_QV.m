clear;
clc;
format long
tic;

%% paramter
myseed = 1;
rng(myseed)

L = 1000;
dt = 1;
Tmax = 1000;
T = 0:dt:Tmax;
nT = length(T);
len = 2^L;
delta0 = -2;
delta = 1;
dk = 2*pi/L;

kx = 2*pi/3-0.1;
ky = 2*pi/(3*sqrt(3));
k = [kx-dk;ky-dk];
kpx = [kx;ky-dk];
kpy = [kx-dk;ky];
d1 = [1;sqrt(3)]/2;
d2 = [1;-sqrt(3)]/2;
d3 = [-1;0];
fk = exp(d1'*k) + exp(d2'*k) + exp(d3'*k);
fkpx = exp(d1'*kpx) + exp(d2'*kpx) + exp(d3'*kpx);
fkpy = exp(d1'*kpy) + exp(d2'*kpy) + exp(d3'*kpy);

%% initial state
Hk0 = [delta0 fk;
    conj(fk) -delta0];
[V,D] = eig(Hk0);
phik0 = V(:,1);
phik = phik0;

Hkpx0 = [delta0 fkpx;
    conj(fkpx) -delta0];
[Vpx,Dpx] = eig(Hkpx0);
phikpx0 = V(:,1);
phikpx = phikpx0;

Hkpy0 = [delta0 fkpy;
    conj(fkpy) -delta0];
[Vpy,Dpy] = eig(Hkpy0);
phikpy0 = Vpy(:,1);
phikpy = phikpy0;

gkk = cell(nT,1);
Gkk = zeros(nT,1);
dphikx = (phikpx0 - phik0)/dk;
dphiky = (phikpy0 - phik0)/dk;
gkk{1} = [dphikx'*dphikx-abs(dphikx'*phik)^2,dphikx'*dphiky-(dphikx'*phik)*(phik'*dphiky);
       dphiky'*dphikx-(dphiky'*phik)*(phik'*dphikx),dphiky'*dphiky-abs(dphiky'*phik)^2];
Gkk(1) = det(gkk{1});

%% time evolution
Hk = [delta fk;
    conj(fk) -delta];
expHk = expm(-1i*Hk*dt);
Hkpx = [delta fkpx;
    conj(fkpx) -delta];
expHkpx = expm(-1i*Hkpx*dt);
Hkpy = [delta fkpy;
    conj(fkpy) -delta];
expHkpy = expm(-1i*Hkpy*dt);

for i = 2:nT
    phik = expHk*phik;
    phikpx = expHkpx*phikpx;
    phikpy = expHkpy*phikpy;
    dphikx = (phikpx - phik)/dk;
    dphiky = (phikpy - phik)/dk;
    gkk{i} = [dphikx'*dphikx-abs(dphikx'*phik)^2,dphikx'*dphiky-(dphikx'*phik)*(phik'*dphiky);
             dphiky'*dphikx-(dphiky'*phik)*(phik'*dphikx),dphiky'*dphiky-abs(dphiky'*phik)^2];
    Gkk(i) = det(gkk{i});
end


%% calculate observable

figure;
% set(gcf, 'position', [250 70 1400 900]);
% % subplot(2,1,1)
plot(T,Gkk);

toc;
