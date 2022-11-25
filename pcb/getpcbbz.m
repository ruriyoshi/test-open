function [ok, bz, rpos, zpos, p_ch] = getpcbbz(rawdata, coeff,date)
% データの並べ替え（積分器のチャンネルとデジタイザの対応）
d=reshape([1:32],2,16);%デジタイザの順番
c=[d;d+32];%積分器の順番
c=c(:);
d_ch=[c; c+64]';

[ch, c2d, d2c] = unique(d_ch);
%d2cはデジタイザ順の配列を積分器順に直し,c2dはその逆
clear c d

%% 
% 無視するチャンネル(積分器＆デジタイザ)

%0~50usの中で分散がゼロはショートしている？
ng1=var(rawdata(1:50,:),1);
%x(:,ng1~=0)

%400~600usの中で明らかに巨大な信号はノイズ？
%%sum(x(411:600,:),1)で400~600usの総和
%外れ値検索：'mean'で
%ng2 = isoutlier(sum(rawdata(411:600,:),1),'mean');
ng2 =sum(abs(rawdata(420:490,:)),1)>=10;
ng4 = sum(abs(rawdata(420:490,:)),1) <=0.3;% | sum(abs(rawdata(420:490,:)),1)>=10;
%構成係数の外れ値を検出して積分器のおかしいチャンネルを抜く
ng3= isoutlier(coeff,'median');
%ok= true(1,128)%& ng2==0 & ng3 == 0;
ok= ng1~=0 & ng2==0 & ng3 == 0 & ng4==0;
clear ng1 ng2 ng3 ng4
 
% 校正係数をかけてbzに直す。今回の構成係数はデジタイザごとになっているので変換の必要はない

bz=rawdata.*coeff;
clear rawdata coeff
%% 
% ch位置の読み込みとチャンネルと位置の対応

r_ch=25;
rind=0.055+0.005*[1,5,9,13,17,[19:41],43,47];
%load('r_pcb2019.mat')
[rpos,zpos] = meshgrid(rind(1:r_ch),-0.042:0.021:0.042);
%p_ch(1,:)=[1:125];

for i=1:5
    p_ch(i,:)=[1:r_ch]+r_ch*(i-1);
end
p_ch(p_ch>=63)=p_ch(p_ch>=63)+1;
%p_ch(p_ch==63)=127;
%p_ch(p_ch==64)=p_ch(p_ch==63);

[ch_p, p2c, c2p] = unique(p_ch(:));

rpos=rpos(p2c);
zpos=zpos(p2c);
clear rind
%% 
% 
% 
% 無視するチャンネル(プローブ)
if date<211212 
%ng_ch=[1,2,5,7,9,13,16,21,[22,25]+25,[2,16]+50,[1,4,13,22]+75];
%ng_ch=[2,[22,25]+25,16+50,[1,4,12,22]+75];
%ng_ch=[69];
ng_ch=[1,13,16,22,47,50,53,64-1];%~1213
%ng_ch=[1,13,16,22,47,50,53];%1214~
%ng_ch=[1,4,13,16,22,[22,25]+25,16+50,2+75,[1,2,5,7,9,13,16,21]+100,[83,95,99,107,111,115,123]-1];
ng=true(1,size(ch_p,1));
ng(ng_ch)=false;
%okp=true(1,128)
okp= ng==true;
%okp=okp(c2p);
clear ng_ch
%% 
% チャンネル変換
%c2d(ch(p_ch));
ok=ok(d2c);
ok=ok(ch(ch_p));
ok=ok&okp;
%ok([15,41,44,52,63,68,88,91,92,119,125])=false;
%ok([66,87,107])=true;

elseif date>=211212 
    load('ok211212.mat');
    if  date>=211214
        ok(63)=true;
    end
end
bz=bz(:,d2c);
bz=bz(:,ch(ch_p));
clear d2c ch_p okp c2p d2c d_ch p2c 
%% 
% Bzの極性をそろえる。440usの時点でのプラスマイナスで判定

bz(:,bz(440,:)<=0)=-bz(:,bz(440,:)<=0);


end