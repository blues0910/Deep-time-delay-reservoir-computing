function [xe,ye] = Equilibrium(tstep,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer,Nv)
%EQUILIBRIUM 此处显示有关此函数的摘要
%   此处显示详细说明
Input_Mask=cell(1,NumberOfLayer);
N=fix(sum(delayOfLayer)/tstep);
x0=rand(N,1);
% x0=rand(Nv*NumberOfLayer,1);
y0=zeros(N,1);
% y0=zeros(Nv*NumberOfLayer,1);
j=0;
for i=1:NumberOfLayer
    j=j+fix(delayOfLayer(i)/tstep);
    y0(j)=(x0(j)-x0(j-1))/tstep;
%     j=Nv*i;
%     y0(j)=(x0(j)-x0(j-1))/tstep;
    Input_Mask{i}=0;
end
i=0;
x=0;
y=0;
temp=1;
while temp>=0.0000001
    i=i+1;
    if i==1
        [xx,x,y]=update_reservior_states(x0,y0,zeros(N,1),zeros(N,1),tstep,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer,Input_Mask,Nv);
    else
        x0=xx;
        y0=y;
        [xx,x,y]=update_reservior_states(x,y,zeros(N,1),zeros(N,1),tstep,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer,Input_Mask,Nv);
        temp=sum(abs(xx-x0).^2);
    end
end
xe=zeros(1,NumberOfLayer);
ye=zeros(1,NumberOfLayer);
j=0;
for i=1:NumberOfLayer
    j=j+fix(delayOfLayer(i)/tstep);
    if i==1
        if deltaOfLayer(i)>0
            xe(i)=0;
        else
            xe(i)=x(j);
        end
        ye(i)=y(j);
    else
        if deltaOfLayer(i)>0
            xe(i)=0;
        else
            xe(i)=x(j);
        end
        ye(i)=betaOfLayer(i)/deltaOfLayer(i)*(sin(kappaOfLayer(i)*xe(i-1)+bOfLayer(i)))^2;
    end
end

% j=0;
% for i=1:NumberOfLayer
%     j=j+fix(delayOfLayer(i)/tstep);
%     if abs(x(j))<=0.0001
%         xe(i)=0;
%     else
%         xe(i)=x(j);
%     end
%     ye(i)=y(j);
% end
end

