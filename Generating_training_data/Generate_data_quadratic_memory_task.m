clear
Q=[0.2 0.1 0.1 0.1; 0.1 0.3 0.2 0.1; 0.1 0.2 0.4 0.1; 0.1 0.1 0.1 0.8];
T=10000;
tmpu=normrnd(0,1,[1,T+length(Q)-1]);
y=zeros(1,T);
for i=1:T
    y(i)=tmpu(1+i-1:4+i-1)*Q*tmpu(1+i-1:4+i-1)';
end
u=tmpu(4:end);
save('3_lag_quadratic_memory_task.mat','u','y','Q')