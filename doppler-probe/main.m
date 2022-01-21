function [] = main(inversion_method,snr,draw_phantom,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare)
%逆変換(1~5),SN(dB),ファントム描画,スペクトル描画,結果描画,等高線(1)/三次元(2),分析描画,比較描画

%パラメータを定義
run define/parameter.m

%軸を定義
run define/axis.m

%描画用ファントムdrF、計算用ファントムFを生成
run make/phantom.m

%ファントムを描画
if draw_phantom
    run draw/phantom.m
end

%W(F->P)を作成
run make/funcW.m

%Pを計算、ノイズを付加
run make/dataP.m

% %Pを保存
% save('mat/save_spectra1.mat','P')

% %Pを読み込む
% load('mat/save_spectra1.mat','P')

% %Pの移動平均をとる
% P = movmean(P, 3);

%観測スペクトルを描画
if draw_spectra
    run draw/spectra.m
end

%---------------Invertion theory----------------

%pinvを使ってW^(-1)を求める。
if inversion_method == 1
    run method/pseudo.m
    run method/goto0.m
end

%Tikhonov 0th
if inversion_method == 2
    run method/Tikhonov0.m
    run method/goto0.m
end

%Tikhonov 1st
if inversion_method == 3
    run method/Tikhonov1.m
    run method/goto0.m
end

%Tikhonov 2nd
if inversion_method == 4
    run method/Tikhonov2.m
    run method/goto0.m
end

%非負制約SIRT
if inversion_method == 5
    run method/SIRT.m
end

%minimum Fischer information
if inversion_method == 6
    run method/MFI.m
    run method/goto0.m
end

%結果の表示
if inversion_method ~= false
    drReF = reshape(ReF, [numVx,numVy]);%描画用に整形
    E = drReF - drF;%エラー=再構成結果-ファントム
    
    %再構成結果、誤差を描画
    if draw_result
        run draw/result.m
    end
    
    %再構成結果から得られるスペクトルを計算、ファントムと比較
    if draw_compare
        run draw/compare.m
    end
end

end
