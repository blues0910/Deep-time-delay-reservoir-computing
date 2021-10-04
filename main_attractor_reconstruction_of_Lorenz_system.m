clear
Training_steps=5000;
Predicting_steps=1000;
T=Training_steps+Predicting_steps;
discarded_steps=100;
training_data=load([pwd '\Generating_training_data\Sample_of_Lorenz_system.mat']);
Input_streaming=training_data.u(:,1:Training_steps);
Training_data=training_data.y(:,1:T);
%------------------------------------------------------------------------------------------------------------
NumberOfLayer=5;
delayOfLayer=[80 40 20 10 10];
deltaOfLayer=[0 0.01 0.01 0.01 0.01];
betaOfLayer=[0.68 0.8 0.97 0.83 0.2];
kappaOfLayer=[0.4 0.1 0.1 0.1 0.1];
bOfLayer=[0.2 0.2 1.5 1.28 1.9];
h=0.2;

Nv=fix(delayOfLayer/1);

x = Equilibrium(h,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer,Nv);
%------------------------------------------------------------------------------------------------------------
Input_Mask=cell(1,NumberOfLayer);
for k=1:NumberOfLayer
    if k==1
        Input_Mask{k}=normalized_Mask(1,size(Input_streaming,1),fix(delayOfLayer(k)/h));
    else
        Input_Mask{k}=ones(1,fix(delayOfLayer(k)/h));
    end
end

epsilon=0;
noise=0+sqrt(epsilon).*randn(fix(sum(delayOfLayer)/h),T);
%------------------------------------------------------------------------------------------------------------
xx=zeros(sum(Nv),T);
X=zeros(fix(sum(delayOfLayer)/h),T);
Y=zeros(fix(sum(delayOfLayer)/h),T);

for it=1:Training_steps
    if it==1
        X0=zeros(fix(sum(delayOfLayer)/h),1)+rand;
        Y0=zeros(fix(sum(delayOfLayer)/h),1)+rand;
    else
        X0=X(:,it-1);
        Y0=Y(:,it-1);
    end
    [xx(:,it),X(:,it),Y(:,it)]=update_reservior_states(X0,Y0,Input_streaming(:,it),noise(:,it),h,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer,Input_Mask,Nv);
end

test=discarded_steps+1;
lamda=0.0000001;
W_out=Batch_training(xx(:,test:Training_steps),Training_data(:,test:Training_steps),lamda);
% W_out=OnLine_training(X(:,test:Training_steps),Training_data(:,test:Training_steps),lamda);

for it=Training_steps+1:T
    X0=X(:,it-1);
    Y0=Y(:,it-1);
    Input=W_out'*xx(:,it-1);
    [xx(:,it),X(:,it),Y(:,it)]=update_reservior_states(X0,Y0,Input,noise(:,it),h,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer,Input_Mask,Nv);
end

y_trained=W_out'*xx(:,1:end);

result=CalcPerf(Training_data,y_trained);
nrmse=result.NRMSE;
rmse=result.RMSE;
mse=result.MSE;
fprintf(' %10.6f',nrmse);
fprintf('\n');

tmpT=discarded_steps+1;
figure('name','1')
tl=linspace(0+(tmpT-1)*h,h*T,size(y_trained,2)-tmpT+1);
for i=1:size(Training_data,1)
    s=[num2str(size(Training_data,1)) '1' num2str(i)];
    subplot(s)
    hold on
    plot(tl(1:Training_steps-tmpT+1),y_trained(i,tmpT:Training_steps),'m')
    plot(tl(Training_steps-tmpT+2:end),y_trained(i,Training_steps+1:end),'r')
    plot(tl,Training_data(i,tmpT:end),'b')
    plot(linspace((Training_steps)*h,(Training_steps)*h,100),linspace(min(y_trained(i,tmpT:end)),max(y_trained(i,tmpT:end)),100),'k')
    xlabel('$t$','Interpreter','latex');
    ylabel('Output','Interpreter','latex');
    legend({'Reservoir output','Predictive output','Actual'},'Interpreter','latex')
end
figure('name','2')
subplot(211)
plot3(y_trained(1,tmpT:end),y_trained(2,tmpT:end),y_trained(3,tmpT:end))
subplot(212)
plot3(Training_data(1,tmpT:end),Training_data(2,tmpT:end),Training_data(3,tmpT:end))



