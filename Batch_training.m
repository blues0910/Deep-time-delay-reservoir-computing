function w_output= Batch_training(xx,y_teach,lamda)
%LNEARREGRESS 此处显示有关此函数的摘要
%   此处显示详细说明
%Tikhonov regularisation or ridge regression---------------------------------------------------
w_output=(y_teach*(pinv(xx'*xx+lamda*diag(var(xx)))*xx'))';
end