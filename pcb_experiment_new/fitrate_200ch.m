%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブの合体率計算
%x点での抵抗率、位置、電場電流密度計算
%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_pr, fitrate, xJt, xEt, xeta, xpos,max_psi_left_pos, max_psi_right_pos] = fitrate_200ch(data2D,grid2D,shot,n,trange)

%
frame=numel(trange);

%max_psi:各時間、各列(z)ごとのpsiの最大値を要素とするvector(dim 1 = z方向(列))。size(1,z,t)
%max_psi_r:各時間、各列の最大値のrのindexを要素とするvector。size(1,z,t)

[max_psi,max_psi_r]=max(data2D.psi,[],1);

%disp(size(max_psi));
%disp(size(max_psi_r));

% squeeze:1の次元を削除
max_psi=squeeze(max_psi);

%disp(size(max_psi));

psi_pr=zeros(3,frame);
xJt=zeros(1,frame);
xEt=xJt;
xpos=zeros(2,frame);%X点座標(z,r)
max_psi_left_pos=zeros(2,frame);%磁気軸座標 psi_pr(1,:)
max_psi_right_pos=zeros(2,frame);%磁気軸座標 psi_pr(2,:)
n_array= 1:n;

for i=1:frame
    r_ind=max_psi_r(:,:,i);

%   islocalmax:各行(2)の局所的最大値, MaxnumExtremaで最も大きい最大値指定、0・1要素の行列を返す
%   find:0出ない要素のindexを返す

    max_psi_ind=find(islocalmax(smooth(max_psi(:,i)),'MaxNumExtrema', 2));
    
    max_psi_ind_left=min(max_psi_ind);
    %disp(max_psi_ind_left);
    max_psi_ind_right=max(max_psi_ind);
    %disp(max_psi_ind_right);

%   min(max_psi_ind):複数max_psi点があった場合小さい方のz座標(左側)
%   max(max_psi_ind):複数max_psi点があった場合大きい方のz座標(右側)

    if numel(max_psi(min(max_psi_ind),i))==0
    psi_pr(1,i)=max_psi(1,i);
    %psi_pr(1,i)=NaN;
    %磁気軸座標の保存
    %disp(i)
    max_psi_left_pos(1,i) = grid2D.zq(1,1);
    max_psi_left_pos(2,i) = grid2D.rq(max_psi_r(1,1,i),1);
    else
    psi_pr(1,i)=max_psi(min(max_psi_ind),i);
    %磁気軸座標の保存
    %disp(i)
    max_psi_left_pos(1,i) = grid2D.zq(1,max_psi_ind_left);
    max_psi_left_pos(2,i) = grid2D.rq(max_psi_r(1,max_psi_ind_left,i),1);
    end
  
    disp(strcat('i=',num2str(i)));
    if numel(max_psi(max(max_psi_ind),i))==0
    psi_pr(2,i)=max_psi(n,i);
    %psi_pr(2,i)=NaN;

    %磁気軸座標の保存
    max_psi_right_pos(1,i) = grid2D.zq(1,n);
    max_psi_right_pos(2,i) = grid2D.rq(max_psi_r(1,n,i),1);
    else
    psi_pr(2,i)=max_psi(max(max_psi_ind),i);
    %磁気軸座標の保存
    max_psi_right_pos(1,i) = grid2D.zq(1,max_psi_ind_right);
    max_psi_right_pos(2,i) = grid2D.rq(max_psi_r(1,max_psi_ind_right,i),1);
    end
    


    if numel(find(islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1)))==0
        psi_pr(3,i)=NaN;
        xJt(1,i)=NaN;
        xEt(1,i)=NaN;
        xpos(:,i)=NaN;
    else
        min_psi_ind=islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1);
        xr=r_ind(min_psi_ind);
       % disp(i)
        xz=n_array(min_psi_ind);
        if xr==1 || xr==n %r両端の場合は検知しない
            psi_pr(3,i)=NaN;
            xJt(1,i)=NaN;
            xEt(1,i)=NaN;
            xpos(:,i)=NaN;
        else
        psi_pr(3,i)=max_psi(min_psi_ind,i);
%         xJt(1,i)=data2D.Jt(xr,min_psi_ind,i);
%         xEt(1,i)=data2D.Et(xr,min_psi_ind,i);
%         xJt(1,i)=min(data2D.Jt(xr-1:xr+1,xz-1:xz+1,i),[],'all');
%         xEt(1,i)=min(data2D.Et(xr-1:xr+1,xz-1:xz+1,i),[],'all');
        %disp(xz);
        %disp(xr);
        xJt(1,i)=mean(data2D.Jt(xr-1:xr+1,xz-1:xz+1,i),'all');       
        xEt(1,i)=mean(data2D.Et(xr-1:xr+1,xz-1:xz+1,i),'all');
        xpos(1,i)=grid2D.zq(1,min_psi_ind);
        xpos(2,i)=grid2D.rq(xr,1);
        
        end
    end



% %左右のprivate fluxの中心
% x1=grid2D.zq(1,min(max_psi_ind));
% y1=grid2D.rq(max_psi_r(1,min(max_psi_ind),i),1);
% x2=grid2D.zq(1,max(max_psi_ind));
% y2=grid2D.rq(max_psi_r(1,max(max_psi_ind),i),1);


%     if max_psi(1,i)==max(max_psi(1:end/2,i))
%         psi_pr(1,i)=NaN; %max_psi(1,i);
%     end    
%     if max_psi(end,i)==max(max_psi(end/2:end,i))
%         psi_pr(2,i)=NaN; %max_psi(end,i);
%     end
end

fitrate=psi_pr(3,:)./min(psi_pr(1:2,:),[],1);

%trange=400:600のときt=430のオフセット引いて合体率を0初めにする
%fitrate=(psi_pr(3,:)-mean(psi_pr(3,1:10),"omitnan"))./min(psi_pr(1:2,:),[],1);
%fitrate = fitrate - fitrate(:,30);
xeta=xEt(1,:)./xJt(1,:);

% %23012763 %456us
% psi_pr(:,20)=NaN;
% fitrate(20)=NaN;
% xJt([19 25])=NaN;
% xEt([1 8 12 19 25])=NaN;
% xeta([1 8 12 19 25 35])=NaN;
% %xeta([5 20 22])=NaN;
% %23011114 
% psi_pr(:,1:3)=NaN;
% fitrate(1:3)=NaN;




end