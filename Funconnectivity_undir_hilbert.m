% this code used the band-pass filter & hilbert transform to extract the
% information of amplitude and phase ,and then computede the functional
% connectivity matrix,such as COH,IMCOH,PLI,PLV,WPLI

close all; clear all; clc;
% %%% this part is not neceassary %%%
% restoredefaultpath; % restore the default of matlab
% toolbox_path = '';
% addpath(genpath(ft_path));

% load the eeg data
Input_datapath = 'E:\Wxj硕博课题\Project01\预处理数据\健康对照组\rest\open\Final3min\';
Input_data = dir(fullfile(Input_datapath,'*.set')); 
FileNames = sort_nat(cellstr({Input_data.name}.'));

% save functional connectivity matrix of each subject
COHsavepath = 'E:\Wxj硕博课题\Project01\三组比较结果\rest\FC_undir_hilbert\COH\open\HC\';
IMCOHsavepath = 'E:\Wxj硕博课题\Project01\三组比较结果\rest\FC_undir_hilbert\IMCOH\open\HC\';
PLVsavepath =  'E:\Wxj硕博课题\Project01\三组比较结果\rest\FC_undir_hilbert\PLV\open\HC\';
PLIsavepath = 'E:\Wxj硕博课题\Project01\三组比较结果\rest\FC_undir_hilbert\PLI\open\HC\';
WPLIsavepath = 'E:\Wxj硕博课题\Project01\三组比较结果\rest\FC_undir_hilbert\WPLI\open\HC\';

%%%%%%%%%%%%%%%%%%%%% compute COH & IMCOH (coherence) %%%%%%%%%%%%%%%%%%%%%
for subj = 1:numel(FileNames)  
    % 导入eeg数据
    EEG = pop_loadset('filename',FileNames{subj},'filepath',Input_datapath);
    % 导入数据重新分段
    EEG = eeg_regepochs(EEG, 'recurrence', 4, 'limits',[0 4], 'rmbase',NaN,'eventtype','o');
    % FFT参数设置
    N = EEG.pnts;
    SampleRate = EEG.srate;
    NFFT = 2^nextpow2(N);
    Freq = SampleRate/2*linspace(0,1,NFFT/2+1);
    for chan = 1:size(EEG.data,1)
        for epochs = 1:size(EEG.data,3)
            %ffts(:,chan,epochs) = fft(hanning(N).*squeeze(EEG.data(chan,:,epochs))',NFFT);% 加窗
            ffts(:,chan,epochs) = fft(squeeze(EEG.data(chan,:,epochs))',NFFT); % 不加窗
        end
    end
    % (1)静息态数据的加窗不同于短时傅里叶变换的加窗
    %  加窗的目的只是为了防止计算功率的时候数据两端的泄露（加不加都行）
    % （2）计划求得的频率点不等同于NFFT的点数，NFFT主要是为了提高频率分辨率
    coh = []; icoh = [];
    for x = 1:size(EEG.data,1)
        for y = 1:size(EEG.data,1)
            fx = squeeze(ffts(:,x,:));
            Pxx = fx.*conj(fx)/N;  % x电极的power
            MeanPx = mean(Pxx,2);  %%%% 平均所有epoch %%%%
            fy = squeeze(ffts(:,y,:));
            Pyy = fy.*conj(fy)/N;  % y电极的power
            MeanPy = mean(Pyy,2); 
            Pxy = fx.*conj(fy)/N;  % Sxy，上述两个电极的交叉谱
            MeanPxy = mean(Pxy,2); 
            C = (abs(MeanPxy).^2)./(MeanPx.*MeanPy); % COH
            IC = imag(MeanPxy) ./ sqrt(MeanPx .* MeanPy); % ICOH
            coh(:,x,y) = C;          % coherence，频率*电极*电极
            icoh(:,x,y) = IC;
        end
    end
    % 清除不需要的频率点
    coh = coh(1:NFFT/2 + 1,:,:);
    icoh = icoh(1:NFFT/2 + 1,:,:);
    % 计算某个频段的平均coh
    frequencyBands = [1 4; 4 8; 8 13; 13 30; 30 50;50 100];
    for freqi = 1:1:size(frequencyBands,1) 
        freqband_indx = frequencyBands(freqi,:);
        Indx = find(freqband_indx(1)<= Freq & Freq <=freqband_indx(2));
        Coh = squeeze(mean(coh(Indx,:,:),1));
        ICoh = squeeze(mean(icoh(Indx,:,:),1));
        COH(:,:,freqi) = Coh;
        IMCOH(:,:,freqi) = ICoh;
    end
    % generate a file name and use the function 'sprintf' to add subj
    [~, file_name, file_extension] = fileparts(EEG.filename);
    filename = sprintf('%s.mat', file_name);
    % using fullfile to consolidate the path and name
    COHpath = fullfile(COHsavepath, filename);
    IMCOHpath = fullfile(IMCOHsavepath, filename);
    %  save 
    save(COHpath, 'COH');
    save(IMCOHpath, 'IMCOH');
    
end
  

%%%%%%%%%%%%%%%% compute PLI & PLV (phase synchronization) %%%%%%%%%%%%%%%%   
for subj = 1:numel(FileNames)
    % 导入eeg数据
    EEG = pop_loadset('filename',FileNames{subj},'filepath',Input_datapath);
    % 定义要分析的频率范围
    frequencyBands = [1 4; 4 8; 8 13; 13 30; 30 50;50 100];
    for freqi = 1:1:size(frequencyBands,1)
        % 数据滤波：对载入的数据进行带通滤波
        % 备注：若为已分段数据，则需调整数据的维度（一般：三维-二维）
        eeg_filtered = eegfilt(EEG.data,EEG.srate,frequencyBands(freqi,1),frequencyBands(freqi,2)); 
        for channels = 1:size(EEG.data,1)
            %逐个电极进行Hilbert变换，并提取相位
            band_phase(channels,:) = angle(hilbert(eeg_filtered(channels,:)));
        end
        % 数据重新分段:定义要分段的大小和数目，电极*样本点*分段
        band_phase = reshape(band_phase,[60 1000 45]); 
        for x = 1:size(band_phase,1)
            for y = 1:size(band_phase,1)
                for epochs = 1:size(band_phase,3)
                    x_phase = squeeze(band_phase(x,:,epochs)); % 提取电极对中第一个电极在某个分段的相位
                    y_phase = squeeze(band_phase(y,:,epochs)); % 提取电极对中第二个电极在某个分段的相位
                    rp = x_phase - y_phase;                    % 计算两个电极在某个分段的相位差
                    % PLV
                    sub_plv(x,y,epochs) = abs(sum(exp(1i*rp))/length(rp));    % 计算被试某个电极对在某个分段的PLV
                    % PLI
                    sub_pli(x,y,epochs) = abs(mean(sign((abs(rp)- pi).*rp))); % 计算被试某个电极对在某个分段的PLI
                end
             end
        end
        pli = mean(sub_pli,3);      % 被试各个分段的PLI计算平均值；pli变量维度是电极*电极
        plv = mean(sub_plv,3);      % 被试各个分段的PLV计算平均值；plv变量维度是电极*电极
        % 不同频率FC矩阵
        PLI(:,:,freqi) = pli;
        PLV(:,:,freqi) = plv;
    end
    % generate a file name and use the function 'sprintf' to add subj
    [~, file_name, file_extension] = fileparts(EEG.filename);
    filename = sprintf('%s.mat', file_name);
    % using fullfile to consolidate the path and name
    PLIpath = fullfile(PLIsavepath, filename);
    PLVpath = fullfile(PLVsavepath, filename);
    %  save 
    save(PLIpath, 'PLI');
    save(PLVpath, 'PLV');
end


%%%%%%%%%%%%%%%%%% compute WPLI (phase synchronization) %%%%%%%%%%%%%%%%%%%
for subj = 1:numel(FileNames)
    % 导入eeg数据
    EEG = pop_loadset('filename',FileNames{subj},'filepath',Input_datapath);
    % 定义要分析的频率范围
    frequencyBands = [1 4; 4 8; 8 13; 13 30; 30 50;50 100];
    for freqi = 1:1:size(frequencyBands,1)
        % 数据滤波：对载入的数据进行带通滤波
        % 备注：若为已分段数据，则需调整数据的维度（一般：三维-二维）
        eeg_filtered = eegfilt(EEG.data,EEG.srate,frequencyBands(freqi,1),frequencyBands(freqi,2)); 
        for channels = 1:size(EEG.data,1)
            % 对每个通道的信号进行hilbert变换，得到解析信号（a + b*i）计算wPLI不需要提取相位
            band_hilbert(channels,:) = hilbert(eeg_filtered(channels,:)); 
        end 
        % 数据重新分段:定义要分段的大小和数目，电极*样本点*分段
        band_hilbert = reshape(band_hilbert,[60 1000 45]);
        for x = 1:size(band_hilbert,1)
            for y = 1:size(band_hilbert,1)
                for epochs = 1:size(band_hilbert,3)
                    x_hilbert = band_hilbert(x,:,epochs);
                    y_hilbert = band_hilbert(y,:,epochs);
                    crossspec = x_hilbert.* conj(y_hilbert); % 交叉谱
                    crossspec_imag = imag(crossspec); % 交叉谱的虚部
                    sub_wpli(x,y,epochs) = abs(mean(crossspec_imag))/mean(abs(crossspec_imag)); 
                    % 计算某个被试某个电极对在某个分段的wPLI
                end
             end
        end
        wpli = mean(sub_wpli,3);
        wpli(isnan(wpli)) = 0;
        % 不同频率FC矩阵
        WPLI(:,:,freqi) = wpli;
    end
    % generate a file name and use the function 'sprintf' to add subj
    [~, file_name, file_extension] = fileparts(EEG.filename);
    filename = sprintf('%s.mat', file_name);
    % using fullfile to consolidate the path and name
    WPLIpath = fullfile(WPLIsavepath, filename);
    %  save 
    save(WPLIpath, 'WPLI');
   
end










