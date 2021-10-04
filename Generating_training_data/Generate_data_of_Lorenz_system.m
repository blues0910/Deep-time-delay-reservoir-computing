clear;
sigma=10;
rho=28;
beta=8/3;
t=[0 1000];
theta=0.01;
N=fix((t(2)-t(1))/theta)+1;
x0=[1 1 1];
sol=dde23('LorenzSystem',[],x0,t,[],sigma, rho, beta);
tmpt=linspace(t(1),t(2),N);
tmp=deval(sol,tmpt);
% plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:))
% hold on
% plot3(tmp(1,:),tmp(2,:),tmp(3,:))
u=tmp(:,499+1:end-1);
y=tmp(:,499+2:end);

save('Sample_of_Lorenz_system.mat');

figure(1)
subplot(311)
plot(y(1,:))
subplot(312)
plot(y(2,:))
subplot(313)
plot(y(3,:))
figure(2)
xp=return_map(u(3,:));
plot(xp(1,:),xp(2,:),'.');
