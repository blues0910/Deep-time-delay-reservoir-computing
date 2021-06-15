clear
Training_steps=10000;
Predicting_steps=1000;
T=Training_steps+Predicting_steps;
discarded_steps=500;
training_data=load([pwd '\Generating_training_data\Sample_of_Lorenz_system.mat']);
Input_streaming=training_data.u(1,1:Training_steps);
Training_data=training_data.y(:,1:T);
%------------------------------------------------------------------------------------------------------------
NumberOfLayer=2;
delayOfLayer=[200 400];
deltaOfLayer=[0 0.01];
betaOfLayer=[0.68 0.8];
kappaOfLayer=[4.0 1];
bOfLayer=[0.2 0.2];
h=0.2;
%------------------------------------------------------------------------------------------------------------
var=0.1;
Input_Mask=-var+(var-(-var)).*rand(fix(delayOfLayer(1)/h),size(Input_streaming,1));
J=Prepocessing_of_input(Input_streaming,Input_Mask);

epsilon=0.000001;
noise=normrnd(0,epsilon,[fix(sum(delayOfLayer)/h),T]);
%------------------------------------------------------------------------------------------------------------
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
    X(:,it)=update_reservior_states(X0,Y0,J{it},noise(:,it),h,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer);
    if mod(it,100)==0
        disp(it)
    end
end

test=discarded_steps+1;
lamda=0.00001;
W_out=Batch_training(X(:,test:Training_steps),Training_data(:,test:Training_steps),lamda);
Output_streaming=W_out'*X(:,1:end);

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



