%%%デジタイザを用いた積分器RC較正データの取得

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

dtacq_num=39;
shot=9;
int_ch=1;
fg_ch=128;

x=getMDSdata(dtacq_num,shot,0);
int=x(:,int_ch);
fg=x(:,fg_ch);
int=int-mean(int);
fg=fg-mean(fg);
RC=rms(int)./rms(fg);

t=1:1000;
figure
plot(t,int,'r')
hold on
plot(t,fg,'b')
legend('integrator','FG')
ylim([-1 1])

