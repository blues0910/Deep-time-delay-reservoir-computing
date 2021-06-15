function J = Prepocessing_of_input(I,M)
%MASKINPUT 此处显示有关此函数的摘要
%   此处显示详细说明
% M is a N*L matrix, 
% N=fix(lag/tstep); the numbers of virtual nodes
% L: the dimension of input signals
% J: time-multiplexed input stream J=M*u, M is a mask vector or matrix
[~,j]=size(I);
J=cell(1,j);
for i=1:j
    J{i}=M*I(:,i);
end
end

