
%% example 1: find the statistics summarised in table 5.1 .

[returnData] = func_estimate_statistics();  % please note that this function will need a very long time to run for 257 days. 




%% example 2: calculate the increase probability given the parameters of order flows
% i.e. the function to calculate the conditional probability based onthe equation (5.28), which is the proposed
% analytic solution of the chapter 5.

x = 200;  % let's say the bid order size on top level is 200;
y = 100;  % the ask order size on top level is 100;


paras = [1.0 1.2 1.0 1.2];  % paras = [p01 p0-1 p10 p-10], denoted in the chapter 5 of the thesis, defined in equation (5.12).

p_inc = alg_pdf_fftv7(x,y,paras);  % this is the function used to calculate the conditional probability of mid-price increase. 




