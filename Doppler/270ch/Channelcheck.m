close all
% clear all
% load("230712_base.mat","data")

%各PCのパスを定義
% run define_path.m
setenv("NIFS_path","N:\")%NIFSのresultsまでのパスを入れる
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.IDS288ch=[pathname.NIFS,'/Doppler/Andor/270CH'];

%------【input】---------------------------------------------------
date = 230913;%【input】実験日
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')

read_data = true;%【input】データをascファイルから読み込む

plot_ICCD = true;%【input】ICCD画像をプロット

cal_CH = true;%【input】CHごとのスペクトルを取得

cal_LineInt = false;%【input】線積分イオン温度、発光強度分布を計算
plot_LineInt_result = false;%【input】線積分イオン温度、発光強度分布をプロット

cal_2D = true;%【input】アーベル変換して2次元イオン温度、発光強度分布を計算
plot_2D_result = true;%【input】2次元イオン温度、発光強度分布をプロット
plot_profile = true;%【input】イオン温度R分布をプロット

%-----------------------解析オプション【input】----------------------------
plot_CH_spectra = false;%【input】CHごとのスペクトルをプロット(cal_CH = trueが必要)
plot_LineInt_interp = false;%【input】死んだCHの補間Ti,Emをプロット(cal_LineInt = trueが必要)
plot_2D_interp = false;%【input】補間スペクトルをプロット(cal_2D = trueが必要)
plot_2D_spectra = 'off';%【input】('off','all','good','bad')2次元スペクトル分布をプロット(cal_2D = trueが必要)

hw_lambda = 50;%【input】波長切り出し半幅
hw_ch = 5;%【input】CH方向切り出し半幅
hw_fit = 12;%【input】フィッティング波長切り出し半幅(< hw_lambda)
num_r = 30;%【input】r分割数(比例して計算時間が増える)
hw_lambdaA = 40;%【input】lambdaA半幅(< hw_lambda)
%------------------------------------------------------------------

%物理定数
% switch ICCD.line
%     case 'Ar'
%         lambda0 = 480.6;%471.3;656.3;480.6;587.562;486.133;656.3;486.133;468.57;486.133;nm%線スペクトル波長
%         mass = 39.95;%4.;12.01%イオン質量数
%     otherwise
%         warning('Input error in ICCD.line.')%ICCD.lineの入力エラー
%         return;
% end

%------校正値(1つのファイルにまとめたい)---------------------------------
% separation = [0 15 31 45 57 73 89 105 120 136 152 164 180 196 211 226 240 255];%CHをZ方向で切り分けるための値
% resolution = -1.420387E-12*lambda0^3 - 2.156031E-09*lambda0^2 + 1.250038E-06*lambda0 + 3.830769E-03;0.0037714;...
%     %-0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769
% z = importdata("z_negative.txt")*1e-3;%計測視線Z[m]
% p = importdata("r.txt")*1e-3;%計測視線と中心軸の距離P[m]
% edge = 0.33;%Rの最大値[m](2022/7　解析～　ポテンシャルの範囲より仮定)
% calib = importdata("Ar_calibration.0916_remake.txt");%ICCD校正ファイル
%------------------------------------------------------------------

if read_data
    %データ読み込み(GUI)
    dir = [pathname.IDS288ch,'/20',num2str(date)];%実験日ディレクトリ
    [file,path] = uigetfile('*.asc','Select a file',dir);
    if isequal(file,0)
        disp('User selected Cancel');
        return
    else
        disp(['User selected ', fullfile(path,file)]);
        filename = fullfile(path,file);
    end
    % %データ読み込み(手入力)
    % filename = [dir,'/shot3_PF39_TF-2_EF120_gas0.62_delay468_width3.asc'];
    %スペクトルを取得
    data = importdata(filename);
    data = data(:,2:1025);%1列目は分光データではないので削除
end
tic

%校正ファイル読み込み
% ch = calib.data(:,1);
% center = calib.data(:,2);
% smile = calib.data(:,3);
% relative = calib.data(:,5);
% instru = calib.data(:,6);
% Ti_instru_CH = 1.69e8*mass*(2.*resolution*instru*sqrt(2.*log(2.))/lambda0).^2;

%各CHの波長軸を生成
% lambda = zeros(2*hw_lambda+1,numel(ch));
% idx_l0 = zeros(numel(ch),1);
% px = transpose(linspace(1,1024,1024));
% for i=1:numel(ch)
%     lx=(px-1-smile(i))*resolution+lambda0;%+0.13
%     idx_l0(i) = knnsearch(lx,lambda0);%lambdaの中で最もlambda0に近いセル番号を取得
%     lambda(:,i) = lx(idx_l0(i)-hw_lambda:idx_l0(i)+hw_lambda);
% end
%%
%ICCD生画像を描画(目視で確認用)

% ピーク検出
sum_data = sum(data(350:420,:),1);
[pks,locs] = findpeaks(sum_data,MinPeakDistance=10);
center=locs(locs>=739);

if plot_ICCD
    f1=figure(1);
    f1.Position=[0 300 400 550];
    subplot(1,2,1)
    contourf(px,px,data')
    title('ICCD Raw Image')
    xlabel('X　(lambda) [px]')
    ylabel('Y (position) [px]')
    hold on
    for i=1:numel(center)
        hol_X = linspace(350,420);
        hol_Y = center(i)*ones(100,1);
        plot(hol_X,hol_Y,'r')
        hold on
        % ver_X = smile(i)*ones(100,1);
        % ver_Y = linspace(center(i)-hw_ch,center(i)+hw_ch);
        % plot(ver_X,ver_Y,'r','LineWidth',1)
    end
    % ylim([900 1024])
    hold off

% f2=figure(2);
    subplot(1,2,2)
    %ICCD生画像を波長方向に足し合わせたもののフィッティング(yを決定する)
    % size(data)

    plot(sum_data)
    hold on
    plot(locs,pks,"o")
    % offset_sum = min(movmean(sum_data,10));
    % sum_data = sum_data - offset_sum;

    % f = fit(px,sum_data,'gauss5');
    % coeff5=coeffvalues(f);
    % centerL = [coeff5(2),coeff5(5),coeff5(8),coeff5(11),coeff5(14)];
    % centerL = round(sort(centerL));
    % subplot(2,1,2)
    % plot(f,px,sum_data')
    % title('ICCD Integrated Image')
    % xlabel('X [px]')
    % ylabel('Strength [a.u.]')
    xlim([1 1024])
    % ylim([0 inf])
    view([90 -90]) % 軸の方向を変更する
    hold off
end

toc