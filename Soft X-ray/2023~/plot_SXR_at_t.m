function [EE_high,EE_low] = plot_SXR_at_t(grid2D,data2D,date,shot,t,show_xpoint,show_localmax,start,interval,save,SXRfilename,filter,NL)
% plot SXR emission on psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: date, date of experiment
%   integer: shot, number of shot
%   integer: t, time of interest (us)
%   boolean: show_xpoint, option for showing the x-point
%   boolean: show_localmax, option for showing the local maximum point
%   integer: start, start time (us)
%   integer: interval, interval time of the framing camera (us)
%   boolean: save, option for saving the reconstruction result
%   string: SXRfilename, name of the SXR image file
%   boolean: filter, option for applying non-linear mean (NLM) filter
%   boolean: NL, option for using non-linear reconstruction

% 再構成パラメータを定めます．_newと付け加えているのは，既に.matファイルに保存されたデータと区別するためです．
N_projection_new = 80;
N_grid_new = 80;

% 実行結果（行列）を保存するフォルダの確認
% なければ作成＆計算、あれば読み込みsave
% 単一時間の場合は保存しない？残りの処理が面倒
% if filter & NL
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/NLF_NLR/',num2str(date),'/shot',num2str(shot));
% elseif ~filter & NL
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/LF_NLR/',num2str(date),'/shot',num2str(shot));
% elseif filter & ~NL
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/NLF_LR/',num2str(date),'/shot',num2str(shot));
% else
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/LF_LR/',num2str(date),'/shot',num2str(shot));
% end

if filter && NL
    options = 'NLF_NLR';
elseif ~filter && NL
    options = 'LF_NLR';
elseif filter && ~NL
    options = 'NLF_LR';
else
    options = 'LF_LR';
end

% 非線形再構成では計算精度を落としたい→グリッド数を落とす
% 投影数は変えない？非線形フィルタなら落とす？
% グリッド数以下の投影数だと正則化する意味もなさそう
% もっと投影数落としていいのかな→その分ノイズ落とす？
% 光ファイバーの本数にしてるけどそこまでする必要もなさそう
% savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/',num2str(date),'/shot',num2str(shot));
% savefolder = strcat('/Users/shinjirotakeda/Documents/Github/SXR_diagnostics/result_matrix/',num2str(date),'/shot',num2str(shot));
% savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/',options,'/',num2str(date),'/shot',num2str(shot));
% filepath = strcat('/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/parameters/',options,'.mat');

% 過去に同じTIF画像を再構成している場合，再構成の計算を行わず直接結果画像を表示させようとします．
% exist(path,'dir')はディレクトリが存在するとき7,存在しないとき0を返します．
savefolder = strcat('G:/My Drive/X-ray/Data/SXROUT/result_matrix/',options,'/',num2str(date),'/shot',num2str(shot));
if exist(savefolder,'dir') == 0
    clc_flag = true;
else
    clc_flag = false;
end

% 再構成計算に必要なパラメータを計算するなら読み込む、しない場合も時間の範囲に関しては読み込む
filepath = '/Users/yuleo/Documents/GitHub/test-open/Soft X-ray/2023~/parameters.mat';
if clc_flag
    if isfile(filepath)
        load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range','N_projection', 'N_grid');
        if N_projection_new ~= N_projection || N_grid_new ~= N_grid
            disp('Calcurate parameters given in plot_SXR_at_t.m didnt match existing parameters in parameters.mat. Recalcurating...' );
            clc_parameters(N_projection_new,N_grid_new,filepath);
            load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
        end
    else
        disp('No existing parameters found. Calcurating...');
        clc_parameters(N_projection_new,N_grid_new,filepath);
        load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
    end
else
    load(filepath,'range');
end

number = (t-start)/interval+1;
plot_flag = false;

if clc_flag
%     if date <= 210924
%         [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter);
%         return
%     else
%         [VectorImage2,VectorImage1] = get_SXRImage(date,number,SXRfilename,filter);
%     end
% TIF画像から再構成されるべき画像をベクトル画像形式で取り出します．ここでフィルタ処理が行われます．
[VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter);
% 再構成計算を行います．ここでペナルティ関数の選択が行われます．
    EE_high = clc_distribution(M,K,gm2d1,U1,s1,v1,VectorImage1,plot_flag,NL);
    EE_low = clc_distribution(M,K,gm2d2,U2,s2,v2,VectorImage2,plot_flag,NL);
else
    % 過去に一度でも再構成を行っている場合，過去のX線分布画像を直接参照しそのまま返そうとします．
    loadpath_high = strcat(savefolder,'/',num2str(number),'_high.txt');
    loadpath_low = strcat(savefolder,'/',num2str(number),'_low.txt');
    EE_high = readmatrix(loadpath_high);
    EE_low = readmatrix(loadpath_low);
end
% subplot(1,2,1);
% imagesc(EE_high);
% subplot(1,2,2);
% imagesc(EE_low);
f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.4];
plot_save_SXR(grid2D,data2D,range,date,shot,t,EE_high,EE_low,show_localmax,show_xpoint,save,filter,NL);

end