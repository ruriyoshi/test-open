function a = fig_et_v_1214(shot)
%     cd("C:\Users\Moe Akimitsu\Documents\MATLAB\研究\matlabver1")
%     psi=importpsi("H:\makimitsu\2020AAPPS\"+string(shot)+"psi.csv");
%     xout=importpsi("H:\makimitsu\2020AAPPS\xout.csv");
%     yout=importpsi("H:\makimitsu\2020AAPPS\yout.csv");
%     t=importpsi("H:\makimitsu\2020AAPPS\"+string(shot)+"time.csv");
%    v=t*0;
    %csvで1次元配列になっているpsiのデータを3次元配列に変更（x,y,t）
    psi=reshape(psi,[size(xout,1),size(yout,1),size(t,1)]);
    %MATLABの配列に合わせるためpsiの配列の次元を入れ替え（y,z,t）
    psi_0=permute(psi,[2 1 3]);
    %
    %Jtの計算のため補完されている点を削減
    z_probe=xout(mod([1:size(xout,1)],30) == 1);%Zは30点補間しているのでそれを戻す
    r_probe=yout(mod([1:size(yout,1)],10) == 1);%Rは10点補間しているのでそれを戻す
    %psi=psi_0(mod([1:size(yout,1)],10) == 1,mod([1:size(xout,1)],30) == 1,:);
    z_probe=xout;
    r_probe=yout;
    psi=psi_0;
    %ψのgradientでBr,Bz,Etを計算する

    [Br,Bz,Et]= gradient(psi,z_probe,r_probe,1/1e6);
    %このままだと1/2πrが計算されてないので
    B = 2*pi*repmat(r_probe,1,size(z_probe,1),size(t,1));
    Br=Br./B;
    Bz=Bz./B;
    Et=Et./B;
    %Jtの計算を行う（）。rotの計算はMATLABのcurl関数で余裕
    for i=1:size(t,1)
        [Jt(:,:,i)]= -curl(z_probe,r_probe,Bz(:,:,i),Br(:,:,i))./(4*pi*1e-7);
    end
    clear Bz Br B
    %磁気中性点をプロットするため、ψのz方向最小値を各ｒに対して求める。min関数で次元を２にする
    iframe=3;
    [psimid,mid]=min(psi_0(:,1:numel(xout)-1,:),[],2,"linear");
    %psimidはmidplaneでのpsi
    %midはxoutの位置　
    cd("C:\Users\Moe Akimitsu\Documents\MATLAB\研究\presentation")
    %plotwithpsi(Jt,z_probe,r_probe,iframe,xout,yout,psi_0)
    %hold on
    %plot(transpose(xout(mid(:,:,iframe))),yout)
    %磁気中性線上のψ極大値がx点，極小値がo点になる。
    %ifo=islocalmin(psimid,1);%極小を求める　pはプロミネンス
    %onum=squeeze(sum(ifo,1)); %各時間でのo点の個数
    ytime=repmat(yout,1,50);
    ttime=repmat(transpose(t),231,1);
    for i=1:size(t,1)
        [opoint(:,i),p(:,i)]=islocalmin(psimid(:,:,i),1);
        [xpoint(:,i),~]=islocalmax(psimid(:,:,i),1);
        [~,maxxp(:,i)]=max(psimid(:,:,i),[],1);
    end
    onum=squeeze(sum(opoint,1));
    sum(onum);
    for i=1:size(t,1)-1
        if onum(i)>=1 && onum(i+1)>=1
            v(i+1)=min(yout(opoint(:,i+1)))-min(yout(opoint(:,i)));
        end
    end
    op=[ttime(opoint),ytime(opoint)];
    save('s'+string(shot)+'.mat','op')
    Etmid=squeeze(Et(mid));
    
    %index=mid(:,:,iframe);
    %plot(xout(index(opoint(:,iframe))),yout(opoint(:,iframe)),"ro")
    %p(opoint(:,iframe));%opointでのpsi
    %hold off
    %X点、O点、
    %opoint時間発展
    h(1) = figure;
    %figure
    yyaxis left %第一軸（左）
    plot(ttime(opoint),ytime(opoint),"ro",'DisplayName','o point')%o点時間発展
    hold on
    plot(ttime(xpoint),ytime(xpoint),"b+",'DisplayName','x point')%x点時間発展
    %plot(t(v ~= 0),v(v ~= 0)*10,"ko-",'DisplayName','v_o')%x点時間発展
    plot(t,v*10,"ko-",'DisplayName','v_o(\times 10^-5)')%x点時間発展
    ylim([-0.3 0.3])
    xlim([460 485])

    yyaxis right %第二軸（右）
    plot(t,squeeze(max(Et,[],[1 2])),"r-",'DisplayName','E_t max')%Et_max時間発展
    %plot(t,max(Etmid,[],1),"b-")%Et_x
    plot(ttime(xpoint),Etmid(xpoint),"g^",'DisplayName','E_t at x point')
    ylim([0 1000])
    hold off
    legend('Location','northwest')
    
    
    
    print(string(shot)+'Vopxpet','-dpng')
%    [~,mid,~]=ind2sub(size(psi_0),mid);%線形インデックスを添え字に戻す
%     %磁気面時間発展プロット
%     h(2) = figure;
%     for i=1:10
%         iframe=i;
%         subplot(2,5,i)
%         plotwithpsi(Jt,z_probe,r_probe,iframe,xout,yout,psi_0)
%         hold on
% 
%         index=mid(:,:,iframe);
%         plot(transpose(xout(index)),yout)
%         plot(xout(index(opoint(:,iframe))),yout(opoint(:,iframe)),"ro")
%         hold off
%     end
%     h(3) = figure;
%     for i=1:10
%         iframe=i+10;
%         subplot(2,5,i)
%         plotwithpsi(Jt,z_probe,r_probe,iframe,xout,yout,psi_0)
%         hold on
% 
%         index=mid(:,:,iframe);
%         plot(transpose(xout(index)),yout)
%         plot(xout(index(opoint(:,iframe))),yout(opoint(:,iframe)),"ro")
%         hold off
%     end
%     h(4) = figure;
%     for i=1:10
%         iframe=i+20;
%         subplot(2,5,i)
%         plotwithpsi(Jt,z_probe,r_probe,iframe,xout,yout,psi_0)
%         hold on
% 
%         index=mid(:,:,iframe);
%         plot(transpose(xout(index)),yout)
%         plot(xout(index(opoint(:,iframe))),yout(opoint(:,iframe)),"ro")
%         hold off
%     end
%     savefig(h,'s'+string(shot)+'.fig',"compact")
%     close(h)
    a=shot
end
