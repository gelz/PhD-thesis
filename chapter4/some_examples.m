%% this function provides two examples for the Chapter 3. 
% author: Geliang Zhang, Please contact zhanggel@gmail.com for any
% questions. 


% the second example .

% please note that this funciton might need support of Matlab toolboxes,
% such as signal processing toolbox and maybe other toolboxes as well. We
% can not provide support for the matlab toolboxes.

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


%% fig 2 & fig 3
% this part we show an example of how to call the VRPF and RBVRPF method and
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

% this the VRPF method;
[T_out2,std_T_out2] = ARTimeVaryingFourierSeries_pitch_tracker(data_input',T_01);

% this the RBVRPF method;
[T_out3,std_T_out3] = RB_pitch_tracker_modified_fast(data_input',k,T_01,SNR_desired_input); % k is the level of noise here. Can be modified to SNR ratio...

%% Figure 2: draw the results of VRPF
Pitch_Proposed = (T_out2./fs)*1000;
Pitch_p_err1 = (T_out2+std_T_out2)./fs*1000;
Pitch_p_err2 = (T_out2-std_T_out2)./fs*1000;

figure();
figure_FontSize = 11;
t_tmp2 = (5:length(Pitch_reference(1:end)))*64/16000;
hold on; 
h = plot(t_tmp2,Pitch_reference(5:end),'k',t_tmp2,Pitch_Proposed(3:end),'r'); 

hold on;
plot(t_tmp2,Pitch_p_err1(3:end),'*r');
hold on;
plot(t_tmp2,Pitch_p_err2(3:end),'*r');

title('comparison of T0 estimated using the VRPF method with the true T0 value');
xlabel(' t/ ms');
ylabel(' T0 /ms');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(h,'LineWidth',1.5);
legend('True T0 value','VRPF','SD of VRPF','Location','South');

%% Figure 3: draw the results of RBVRPF

Pitch_Proposed = (T_out3./fs)*1000;
Pitch_p_err1 = (T_out3+std_T_out3)./fs*1000;
Pitch_p_err2 = (T_out3-std_T_out3)./fs*1000;

Pitch_Proposed = downsample(Pitch_Proposed,64);
Pitch_p_err1 = downsample(Pitch_p_err1,64);
Pitch_p_err2 = downsample(Pitch_p_err2,64);

figure();
figure_FontSize = 11;
t_tmp2 = (5:length(Pitch_reference(1:end)))*64/16000;
hold on; 
h = plot(t_tmp2,Pitch_reference(5:end),'k',t_tmp2,Pitch_Proposed(6:end),'r'); 

hold on;
plot(t_tmp2,Pitch_p_err1(6:end),'*r');
hold on;
plot(t_tmp2,Pitch_p_err2(6:end),'*r');

title('comparison of T0 estimated using the RBVRPF method with the true T0 value');
xlabel(' t/ ms');
ylabel(' T0 /ms');
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(h,'LineWidth',1.5);
legend('True T0 value','RBVRPF','SD of RBVRPF','Location','South');

