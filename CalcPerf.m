function Result=CalcPerf(Refernce,Test)
% INPUT 
% Refernce M x N
% Test M x N
% Output
% Result-struct
% 1.MSE (Mean Squared Error)
% 2.PSNR (Peak signal-to-noise ratio)
% 3.R Value
% 4.RMSE (Root-mean-square deviation)
% 5.NRMSE (Normalized Root-mean-square deviation)
% 6.MAPE (Mean Absolute Percentage Error)
% Developer Abbas Manthiri S
% Mail Id abbasmanthiribe@gmail.com
% Updated 27-03-2017
% Matlab 2014a
%% geting size and condition checking
[row_R,col_R,dim_R]=size(Refernce);
[row_T,col_T,dim_T]=size(Test);
if row_R~=row_T || col_R~=col_T || dim_R~=dim_T
    error('Input must have same dimentions')
end
%% Common function for matrix
% Mean for Matrix
meanmat=@(a)(mean(mean(a)));
% Sum for Matrix
summat=@(a)(sum(sum(a)));
% Min  for Matrix
minmat=@(a)(min(min(a)));
% Max  for Matix
maxmat=@(a)(max(max(a)));
%% MSE Mean Squared Error
Result.MSE = meanmat((Refernce-Test).^2);
%% PSNR Peak signal-to-noise ratio
range=[1,255];
if max(Refernce(:))>1
    maxI=range(2);
else
    maxI=range(1);
end
Result.PSNR= 10* log10(maxI^2/Result.MSE);
%% R Value
Result.Rvalue=1-abs( summat((Test-Refernce).^2) / summat(Refernce.^2) );
%% RMSE Root-mean-square deviation
Result.RMSE=abs( sqrt( meanmat((Test-Refernce).^2) ) );
%% Normalized RMSE Normalized Root-mean-square deviation
Result.NRMSE=Result.RMSE/(maxmat(Refernce)-minmat(Refernce));
%% MAPE Mean Absolute Percentage Error
Result.Mape=meanmat(abs(Test-Refernce)./Refernce)*100;

