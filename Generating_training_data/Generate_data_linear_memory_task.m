clear
L=[0.2 0.3 0.5];
T=10000;
tmpu=normrnd(0,1,[1,T+length(L)-1]);
y=zeros(1,T);
for i=1:T
    y(i)=L*tmpu(1+i-1:3+i-1)';
end
u=tmpu(3:end);
save('2_lag_linear_memory_task.mat','u','y','L')