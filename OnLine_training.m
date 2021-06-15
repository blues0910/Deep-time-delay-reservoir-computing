function [W_out] = OnLine_training(x,y,alpha)
%TRAINONLINE 此处显示有关此函数的摘要
%   此处显示详细说明
Ph=eye(size(x,1),size(x,1))/alpha;
W_out=zeros(size(x,1),size(y,1));
for it=1:size(x,2)
    [W_out,Ph]=iterations(W_out,Ph,x(:,it),y(:,it));
end
end

function [W,P] = iterations(Wh,Ph,r,f)
e=Wh'*r-f;
tmp11=Ph*r;
tmp12=tmp11*r';
tmp1=tmp12*Ph;
tmp2=1+r'*Ph*r;
P=Ph-tmp1/tmp2;
W=Wh-P*r*e';
end


