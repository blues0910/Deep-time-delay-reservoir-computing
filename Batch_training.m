function w_output= Batch_training(xx,y_teach,lamda)
%LNEARREGRESS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%Tikhonov regularisation or ridge regression---------------------------------------------------
w_output=(y_teach*(pinv(xx'*xx+lamda*diag(var(xx)))*xx'))';
end