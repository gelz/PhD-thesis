%% this function provides two examples for the Chapter 3. 
% author: Geliang Zhang, Please contact zhanggel@gmail.com for any
% questions. 


% the first example 

%% fig1: draw a clean sample of the speech signal. 

clear all; clc;
[d0,fs] = wavread('mic_M05_si1211.wav');
d1 = downsample(d0,3);
fs = fs/3;
begin1 = 65001;
over1  = 68000;
data_cut1 = d1(begin1:over1);

figure();
plot((1:length(data_cut1))./fs*1000,data_cut1);

figure_FontSize=11;
title('sample of speech','fontsize',figure_FontSize);
xlabel(' t/ ms');
ylabel('Amplitude');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);


%% fig 2
% this part we show an example of how to call the partile filter method and
% compare the results with the true result. 

% clear all; clc;
[d0,fs] = audioread('mic_M05_si1211.wav');  % this is the data file; 
d1 = downsample(d0,3);
fs = fs/3;
begin1 = 65001;
over1  = 68000;
data_cut1 = d1(begin1:over1);
[d2,fs] = audioread('lar_M05_si1211.wav');  % this is the file to extract true pitch period as benchmark. 
d3 = downsample(d2,3);
fs = fs/3;
[fx_rapt,tt_rapt]=fxrapt_newzgl(d3(begin1-640:over1+640),fs); % this is the function to produce the true pitch period using clean speech. 
pitch_true_cut = fx_rapt(6:end-5);
Pitch_reference = 1000./pitch_true_cut;


SNR_desired_input = 0;  % we control the level of noise here. in units of dB. 

data_sample = d1;

noise = randn(length(data_sample),1);
prior_SNR = 10*log10((sum(data_sample((begin1-640):(over1+640)).^2)./sum((noise((begin1-640):(over1+640))).^2)));
k = 10^((prior_SNR-SNR_desired_input)/20);


y_ob = data_sample + k*noise;
y = y_ob;
data_sample2 = data_sample;
data_sample = y;

begin1 = 65001;
over1  = 68000;
clean_data = d1(begin1:over1);

d1 = data_sample;
data_cut1 = d1(begin1:over1);

T_01 = round(32000./(pitch_true_cut(1)+pitch_true_cut(2)));

data_input = data_cut1;
% [T_out2,std_T_out2] =
% PSLOA_model_pitch_tracker(data_input',T_01,clean_data); % we can use the
% first three periods of clean data to stabilize the method. But it can
% work without it. 
[T_out2,std_T_out2] = PSLOA_model_pitch_tracker(data_input',T_01,data_input);


Pitch_Proposed = (T_out2./fs)*1000;
Pitch_p_err1 = (T_out2+std_T_out2)./fs*1000;
Pitch_p_err2 = (T_out2-std_T_out2)./fs*1000;

figure();
figure_FontSize = 11;
t_tmp2 = (5:length(Pitch_reference(1:end)))*64/16000;
hold on; 
h = plot(t_tmp2,Pitch_reference(5:end),'k',t_tmp2,Pitch_Proposed(5:end),'r'); 

hold on;
plot(t_tmp2,Pitch_p_err1(5:end),'*r');
hold on;
plot(t_tmp2,Pitch_p_err2(5:end),'*r');

title('comparison of T0 estimated using the proposed particle filter method with the true T0 value');
xlabel(' t/ ms');
ylabel(' T0 /ms');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(h,'LineWidth',1.5);
legend('True T0 value','Particle Filter','SD of Particle Filter','Location','South');

