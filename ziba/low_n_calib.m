M = readmatrix( "C:\Users\uswk0\OneDrive\デスクトップ\data\low_n_calib\BT&-0100.CSV" );

%%%%%%%%%各パラメータ定義%%%%%%%%

%時間空間でのサンプリング間隔
Ts = 1.00E-07;
%Ts秒ごとに1E-03秒間にわたってサンプリングされた時間ベクトルt
t = 0:1.00E-07:1.00E-03;
%周波数空間でのサンプリング間隔
%sampling frequency
Fs = 1/Ts;
%サンプリング数n
len_freq = length(t);

perc = 0.02; % span used in smoothing
windowSize = 5; % window size in filtering


%移動平均フィルターの有理伝達関数 b/a
b = (1/windowSize)*ones(1,windowSize);
a = 1;

%%%%%%% smoothing and filtering %%%%%%%

%smooth：重みづけで外れ値をなめらかにする
% 'rloess' → 回帰において外れ値に小さい重みを割り当て。
% 'perc'   → '0,02':合計データ点数の 2% の範囲を使用

%filtering:平均をとって細かい振動をなめらかな曲線にする
%今回は移動平均フィルター。データに沿って長さ windowSize のウィンドウをスライド、各ウィンドウに含まれるデータの平均を計算
%b/a有理伝達関数

smoothed = smooth(M(:,3),perc,'rloess');
filtered_and_smoothed = filter(b,a,smoothed);

%%%%%%% fft %%%%%%%

%両側スペクトル P2 を計算 → P2 および偶数の信号長len_freqに基づいて、片側スペクトル P1 を計算。
%周波数空間での信号のサンプリングに対応するベクトル f 

y_fft = fft(M(:,3));
P2 = abs(y_fft/len_freq);
P1 = P2(1:(len_freq-1)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(len_freq/2))/len_freq;

%%%%%%% plotting %%%%%%%%
figure
subplot(2,1,1)
hold on
plot(t,filtered_and_smoothed,'LineWidth',1);
scatter(t,M(:,3),1);
xlabel("Time (t)")
ylabel("Signal (V)")
hold off
subplot(2,1,2)
plot(f(1:100),P1(1:100)) 
xlabel("f (Hz)")
ylabel("Amplitude (V)")

%%%%%%%% calculate rms and peak %%%%%%%
rms_raw = rms(M(:,3));
rms_smoothed = rms(filtered_and_smoothed);
peak_smoothed = rms_smoothed * sqrt(2);
disp(strcat('rms_raw = ', num2str(rms_raw)));
disp(strcat('rms_smoothed = ', num2str(rms_smoothed)));
disp(strcat('peak_smoothed = ', num2str(peak_smoothed)));