function rms_mat = low_n_calib(direction)

% rms_mat:(probe本数) × 4　(rms_func_gene, rms_amp, rms_coil_raw, rms_coil_smoothed)の行列を計算する関数
% direction = 'T' or 'R' or 'Z'


probe_num = 8;
rms_mat = zeros(probe_num,4);

for i = 1:8
    disp(strcat('進行確認 i = ', num2str(i)));
    
    %終端抵抗無しファイル → '--', 終端抵抗ありファイル → '&-'
    filename = strcat('B',direction,'--','0',num2str(i),'00');
    %filename = strcat('B',direction,'&-','0',num2str(i),'00');

    M = readmatrix(strcat("C:\Users\uswk0\OneDrive\デスクトップ\道家\卒論\data\low_n_calib\", filename, '.CSV'));
    %オシロデータの最初の30行無視
    M = M(31:10031,:);

    %%%%%%%%% 各パラメータ定義　%%%%%%%%

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
    % 'perc'   → '0.02':合計データ点数の 2% の範囲を使用

    %filtering:平均をとって細かい振動をなめらかな曲線にする
    %今回は移動平均フィルター。データに沿って長さ windowSize のウィンドウをスライド、各ウィンドウに含まれるデータの平均を計算
    %b/a有理伝達関数

    smoothed = smooth(M(:,3),perc,'rloess');
    filtered_and_smoothed = filter(b,a,smoothed);

    %%%%%%% fft %%%%%%%
    %実験や観測データの解析では0<f<∞とするほうが便利なため，両側スペクトルの大きさを2倍にして0<f<∞の範囲で定義した片側スペクトルを使う
    %両側スペクトル(−∞<f<∞) P2 を計算 → P2 および偶数の信号長len_freqに基づいて、片側スペクトル P1(0<f<∞) を計算。
    %周波数空間での信号のサンプリングに対応するベクトル f 

    %ch3の信号のみ抽出(3列目)
    y_fft = fft(M(:,3));
    P2 = abs(y_fft/len_freq);
    P1 = P2(1:(len_freq-1)/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(len_freq/2))/len_freq;

    %{
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
    savefig(strcat("C:\Users\uswk0\OneDrive\デスクトップ\道家\卒論\data\low_n_calib_picture\",filename));
    
    %}

    %%%%%%%% calculate rms and peak %%%%%%%
    rms_func_gene = rms(M(:,1));
    rms_amp = rms(M(:,2));
    rms_coil_raw = rms(M(:,3));
    rms_coil_smoothed = rms(filtered_and_smoothed);
    peak_coil_smoothed = rms_coil_smoothed * sqrt(2);

    %{
    disp(strcat('rms_func_gene = ', num2str(rms_func_gene)));
    disp(strcat('rms_amp = ', num2str(rms_amp)));
    disp(strcat('rms_coil_raw = ', num2str(rms_coil_raw)));
    disp(strcat('rms_coil_smoothed = ', num2str(rms_coil_smoothed)));
    disp(strcat('peak_coil_smoothed = ', num2str(peak_coil_smoothed)));
    %}

    rms_mat(i,:) = [rms_func_gene, rms_amp, rms_coil_raw, rms_coil_smoothed];
end
