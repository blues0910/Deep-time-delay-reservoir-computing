function [X,Y]=update_reservior_states(X0,Y0,Input,z,h,NumberOfLayer,delayOfLayer,deltaOfLayer,betaOfLayer,kappaOfLayer,bOfLayer,Input_Mask)
% 此处显示有关此函数的摘要
% 此处显示详细说明
%---------------------------
%---------------------------
X=zeros(fix(sum(delayOfLayer)/h),1);
Y=zeros(fix(sum(delayOfLayer)/h),1);
for i=1:NumberOfLayer
    N=fix(sum(delayOfLayer(1:i))/h);
    if i==1
        I1=times(Input_Mask{i},Input);
        [X(1:N),Y(1:N)] = layer(delayOfLayer(i),X0(1:N),Y0(N),h,I1,z(1:N),deltaOfLayer(i),betaOfLayer(i),kappaOfLayer(i),bOfLayer(i));
        if NumberOfLayer>1
            I2=pchip(linspace(1,delayOfLayer(i),N),X(1:N),linspace(1,delayOfLayer(i),fix(delayOfLayer(i+1)/h)));
        end
    else
        NN=fix(sum(delayOfLayer(1:i-1))/h);
        I2=times(Input_Mask{i},I2);
        [X(NN+1:N),Y(NN+1:N)] = layer(delayOfLayer(i),X0(NN+1:N),Y0(fix(sum(delayOfLayer(1:i))/h)),h,I2,z(NN+1:N),deltaOfLayer(i),betaOfLayer(i),kappaOfLayer(i),bOfLayer(i));
        if i~=NumberOfLayer
            I2=pchip(linspace(1,delayOfLayer(i),fix(delayOfLayer(i)/h)),X(NN+1:N),linspace(1,delayOfLayer(i),fix(delayOfLayer(i+1)/h)));
%             I2=pchip(X(NN+1:N),linspace(1,delayOfLayer(i),fix(delayOfLayer(i)/h)),linspace(1,delayOfLayer(i),fix(delayOfLayer(i+1)/h)));
        end
    end
end
end

function [x,y] = layer(delay,x0,y0,h,J,z,delta,beta,kappa,b)
N=fix(delay/h);
x=zeros(N,1);
y=zeros(N,1);
for j=1:N 
    if j==1
        [x(j),y(j)]=marunge4(@Fx,@Fy,x0(N),y0,h,x0(1),J(j),z(j),delta,beta,kappa,b);
    else
        [x(j),y(j)]=marunge4(@Fx,@Fy,x(j-1),y(j-1),h,x0(j),J(j),z(j),delta,beta,kappa,b);
    end
end
end
function [x,y]=marunge4(varargin)
FuncHandle1=varargin{1};
FuncHandle2=varargin{2};
x0=varargin{3};
y0=varargin{4};
h=varargin{5};
xh=varargin{6};
J=varargin{7};
z=varargin{8};
delta=varargin{9};
beta=varargin{10};
kappa=varargin{11};
b=varargin{12};

% k1=FuncHandle1(x0,y0,xh,J,z,delta,beta,kappa,b);
% l1=FuncHandle2(x0);
% 
% k2=FuncHandle1(x0+k1/2,y0+l1/2,xh,J,z,delta,beta,kappa,b);
% l2=FuncHandle2(x0+k1/2);
% 
% k3=FuncHandle1(x0+k2/2,y0+l2/2,xh,J,z,delta,beta,kappa,b);
% l3=FuncHandle2(x0+k2/2);
% 
% k4=FuncHandle1(x0+k3,y0+l3,xh,J,z,delta,beta,kappa,b);
% l4=FuncHandle2(x0+k3);
% 
% x=x0+h*(k1+2*k2+2*k3+k4)/6;
% y=y0+h*(l1+2*l2+2*l3+l4)/6;
x=x0+h*FuncHandle1(x0,y0,xh,J,z,delta,beta,kappa,b);
y=y0+h*FuncHandle2(x0);
end
function v = Fx(varargin)
x=varargin{1};
y=varargin{2};
xh=varargin{3};
J=varargin{4};
z=varargin{5};
delta=varargin{6};
beta=varargin{7};
kappa=varargin{8};
b=varargin{9};

v=-x-y*delta+beta*(sin(xh+kappa*J+b+z))^2;
end
function v = Fy(varargin)
y=varargin{1};
v=y;
end

