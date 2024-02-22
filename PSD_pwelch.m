% this code was used to compute the power spectrum density form the
% preprocessing EEG data based on the pwelch.m function. usually, this code
% fitted the resting-state EEG data,included, absolute psd & relative psd

close all; clear all; clc;
% load the data file
filepath = '\';
Input = dir(fullfile(filepath,'*.set')); 
FileName = sort_nat(cellstr({Input.name}.'));  
for subj = 1:1:length(FileName)
    % load EEG data based on eeglab
    EEG = pop_loadset('filename',FileName{subj},'filepath',filepath);
    Fs = EEG.srate;       % data resample
    L = EEG.pnts;         % data length
    NFFT = 2^nextpow2(L); %计算用于FFT的下一个最接近2的幂次方的数值，以确保FFT计算的效率
    for ii = 1:size(EEG.data,1)
        for jj=1:size(EEG.data,3) % the number of EEG data epochs
            y = squeeze(EEG.data(ii,:,jj));
           [pxx2, freq] = pwelch(y, hanning(512), 256, NFFT, Fs); 
           % 这部分的参数设置根据自己数据的情况决定，一般窗选择汉宁窗，窗长为数据长度的一半，窗的重叠为窗长的一半
           psd2(ii,:,jj) = 10*log10(pxx2);
        end
    end
    psd_avg2 = squeeze(mean(psd2,3)); % average all epochs
    Psd_grp_HCopen(:,:,subj) = psd_avg2;
end

% plot the psd figure
figure
plot(freq,mean(Psd_grp_HCopen1,3)); 

% compute the psd of different frequencies
freqBands = [1 4; 4 8; 8 13; 13 30; 30 50;50 100];
% 频段的选择可参考之前的文献
delta_HCopen = squeeze(mean(Psd_grp_HCopen(:,find(freqBands(1,1)<freq & freq<freqBands(1,2)),:),2));
theta_HCopen = squeeze(mean(Psd_grp_HCopen(:,find(freqBands(2,1)<freq & freq<freqBands(2,2)),:),2));
alpha_HCopen = squeeze(mean(Psd_grp_HCopen(:,find(freqBands(3,1)<freq & freq<freqBands(3,2)),:),2));
beta_HCopen = squeeze(mean(Psd_grp_HCopen(:,find(freqBands(4,1)<freq & freq<freqBands(4,2)),:),2));
lgamma_HCopen = squeeze(mean(Psd_grp_HCopen(:,find(freqBands(5,1)<freq & freq<freqBands(5,2)),:),2));
hgamma_HCopen = squeeze(mean(Psd_grp_HCopen(:,find(freqBands(6,1)<freq & freq<freqBands(6,2)),:),2));
gamma_HCopen = squeeze(mean(Psd_grp_HCopen(:,find(freqBands(5,1)<freq & freq<freqBands(6,2)),:),2));
