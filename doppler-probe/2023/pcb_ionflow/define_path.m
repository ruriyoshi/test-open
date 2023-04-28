%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
setenv("NIFS_path","/Volumes/experiment/results")
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
setenv('rsOnedrive','/Users/rsomeya/Library/CloudStorage/OneDrive-TheUniversityofTokyo/lab')
setenv("rsGdrive","/Users/rsomeya/Library/CloudStorage/GoogleDrive-rsomeya2016@g.ecc.u-tokyo.ac.jp/マイドライブ/lab")
pathname.fig=[getenv('rsOnedrive') '/figure'];%figure保存先
pathname.mat=[getenv('rsOnedrive') '/mat'];%mat保存先
pathname.rawdata=[pathname.mat,'/pcb_raw'];%dtacqのrawdata.matの保管場所
pathname.processeddata=[pathname.mat,'/pcb_processed'];%磁場データmatの保管場所
pathname.flowdata=[pathname.mat,'/ionflow'];%流速データmatの保管場所
pathname.vdistdata=[pathname.mat,'/ionvdist'];%速度分布データmatの保管場所
% pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
% pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
% pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
% pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
% pathname.save=getenv('output');%outputデータ保存先