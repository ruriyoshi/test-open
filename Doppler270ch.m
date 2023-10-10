close all
% clear all
% load("230712_base.mat","data")

%各PCのパスを定義
run define_path.m
%setenv("NIFS_path","N:\")%NIFSのresultsまでのパスを入れる
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.IDS270ch=[pathname.NIFS,'/Doppler/Andor/270CH'];

%------【input】---------------------------------------------------
date = 230920;%【input】実験日
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')

read_data = true;%【input】データをascファイルから読み込む

plot_ICCD = true;%【input】ICCD画像をプロット

cal_CH = true;%【input】CHごとのスペクトルを取得

cal_LineInt = true;%【input】線積分イオン温度、発光強度分布を計算
plot_LineInt_result = true;%【input】線積分イオン温度、発光強度分布をプロット

cal_2D = true;%【input】アーベル変換して2次元イオン温度、発光強度分布を計算
plot_2D_result = true;%【input】2次元イオン温度、発光強度分布をプロット
plot_profile = true;%【input】イオン温度R分布をプロット

%-----------------------解析オプション【input】----------------------------
plot_CH_spectra = true;%【input】CHごとのスペクトルをプロット(cal_CH = trueが必要)
plot_LineInt_interp = false;%【input】死んだCHの補間Ti,Emをプロット(cal_LineInt = trueが必要)
plot_2D_interp = false;%【input】補間スペクトルをプロット(cal_2D = trueが必要)
plot_2D_spectra = 'all';%【input】('off','all','good','bad')2次元スペクトル分布をプロット(cal_2D = trueが必要)

hw_lambda = 50;%【input】波長切り出し半幅
hw_ch = 4;%5;%【input】CH方向切り出し半幅
hw_fit = 12;%【input】フィッティング波長切り出し半幅(< hw_lambda)
num_r = 30;%【input】r分割数(比例して計算時間が増える)
hw_lambdaA = 40;%【input】lambdaA半幅(< hw_lambda)
%------------------------------------------------------------------

%物理定数
switch ICCD.line
    case 'Ar'
        lambda0 = 480.6;%471.3;656.3;480.6;587.562;486.133;656.3;486.133;468.57;486.133;nm%線スペクトル波長
        mass = 39.95;%4.;12.01%イオン質量数
    otherwise
        warning('Input error in ICCD.line.')%ICCD.lineの入力エラー
        return;
end

%------校正値(1つのファイルにまとめたい)---------------------------------
%separation = [0 15 31 45 57 73 89 105 120 136 152 164 180 196 211 226 240 255];%CHをZ方向で切り分けるための値%todo
resolution = -1.420387E-12*lambda0^3 - 2.156031E-09*lambda0^2 + 1.250038E-06*lambda0 + 3.830769E-03;0.0037714;...
    %-0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769
%z = importdata("z_negative.txt")*1e-3;%計測視線Z[m]
%p = importdata("r.txt")*1e-3;%計測視線と中心軸の距離P[m]
edge = 0.36;%0.33;%Rの最大値[m](2022/7　解析～　ポテンシャルの範囲より仮定)
%calib = importdata("Ar_calibration.0916_remake.txt");%ICCD校正ファイル

%較正係数のバージョンを日付で判別
calib_filename='Doppler270ch_Ar_dial4820.xlsx';
sheets = sheetnames(calib_filename);
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));
C = readmatrix(calib_filename,'Sheet',num2str(sheet_date));

z=unique(C(:,7))*1e-3;%計測点のZ座標[m]
%------------------------------------------------------------------
%%
%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','doCheck:'};
dlgtitle = 'Input';
dims = [1 35];
if exist('date','var') && exist('IDXlist','var') && exist('doCheck','var')
    definput = {num2str(date),num2str(IDXlist+1),num2str(doCheck)};
else
    definput = {'','',''};
end
% definput = {'','',''};
% definput = {num2str(date),num2str(IDXlist),num2str(doCheck)};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
doCheck = logical(str2num(cell2mat(answer(3))));
if read_data
    %データ読み込み(GUI)
    dir = [pathname.IDS270ch,'/',num2str(date)];%実験日ディレクトリ
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
%%
screensize=get(0,'screensize');
fig_position=[0,50,screensize(3)-50,screensize(4)-150];

%校正ファイル読み込み
is_alive=find(C(:,9));%生きているチャンネルのインデックスのみ抜き出す
%以下生きているチャンネルしか扱わない
ch=C(is_alive,1);
center=round(C(is_alive,2),0);%calibファイルが小数なので，四捨五入
smile=round(C(is_alive,3))+35;%4175に合わせてsmileを作ったので，手動で定数を足して合わせている
relative=C(is_alive,4);
instru=C(is_alive,5);
p=C(is_alive,8)*1e-3;
separation=cat(1,0,find(ischange(C(is_alive,6)))-1);%そのzで最後のインデックス(notチャンネル番号)

Ti_instru_CH = 1.69e8*mass*(2.*resolution*instru*sqrt(2.*log(2.))/lambda0).^2;

%各CHの波長軸を生成
lambda = zeros(2*hw_lambda+1,numel(ch));
idx_l0 = zeros(numel(ch),1);
px = transpose(linspace(1,1024,1024));
for i=1:numel(ch)
    lx=(px-1-smile(i))*resolution+lambda0;%+0.13
    idx_l0(i) = knnsearch(lx,lambda0);%lambdaの中で最もlambda0に近いセル番号を取得
    lambda(:,i) = lx(idx_l0(i)-hw_lambda:idx_l0(i)+hw_lambda);
end

%ICCD生画像を描画(目視で確認用)
if plot_ICCD
    %figureの設定
    fig_num=1;
    f=figure(fig_num);
    f.Position=fig_position;

    %1枚目(左)
    subplot(1,2,1)
    contourf(px,px,data',50,'LineStyle','none')%変える
    colormap turbo
    colorbar
    
    title('ICCD Raw Image')
    xlabel('X　(lambda) [px]')
    ylabel('Y (position) [px]')
    hold on
    for i=1:numel(ch)
        hol_X = linspace(smile(i)-hw_lambda,smile(i)+hw_lambda);
        hol_Y = center(i)*ones(100,1);
        plot(hol_X,hol_Y,'r')
        hold on
        ver_X = smile(i)*ones(100,1);
        ver_Y = linspace(center(i)-hw_ch,center(i)+hw_ch);
        plot(ver_X,ver_Y,'r','LineWidth',1)
    end
    % ylim([900 1024])
    hold off
    %ICCD生画像をCH方向に足し合わせたもののフィッティング(目視で確認用)
    sum_data = sum(data,2);
    offset_sum = min(movmean(sum_data,10));
    sum_data = sum_data - offset_sum;
    f = fit(px,sum_data,'gauss5');
    coeff5=coeffvalues(f);
    % centerL = [coeff5(2),coeff5(5),coeff5(8),coeff5(11),coeff5(14)];
    % centerL = round(sort(centerL));
    subplot(1,2,2)
    plot(f,px,sum_data')
    title('ICCD Integrated Image')
    xlabel('X [px]')
    ylabel('Strength [a.u.]')
    xlim([1 1024])
    ylim([0 inf])
    hold off
end
%%
%-------CHごとの線積分温度、線積分発光強度を計算--------
if cal_CH
    hw_plot_px = hw_lambda;%プロット範囲半幅
    col_subp1 = 2;%サブプロット行数
    raw_subp1 = 2;%サブプロット列数
    n_subp1 = col_subp1*raw_subp1;
    passive_Ti = zeros(numel(ch),1);%CH温度[eV]
    passive_Timax = zeros(numel(ch),1);%CH温度95%信頼区間上限[eV]
    passive_Timin = zeros(numel(ch),1);%CH温度95%信頼区間下限[eV]
    passive_Em = zeros(numel(ch),1);%CH発光強度[a.u.]
    spectra = zeros(2*hw_lambda+1,numel(ch));%CHスペクトル
    for i=1:numel(ch)
        spectra(:,i) = sum(data(idx_l0(i)-hw_lambda:idx_l0(i)+hw_lambda,center(i)-hw_ch:center(i)+hw_ch),2)/relative(i);
        offset = min(movmean(spectra(:,i),20));
        spectra(:,i) = spectra(:,i) - offset;%オフセットを引く
        f = fit(lambda(hw_lambda+1-hw_fit:hw_lambda+1+hw_fit,i),spectra(hw_lambda+1-hw_fit:hw_lambda+1+hw_fit,i),'gauss1');
        coeff = coeffvalues(f);%フィッティング係数
        % coeff(3)/resolution/sqrt(2);%IDLのcoeff[2]に等しい(確認用)
        % Sigma(i) = coeff(3)/sqrt(2);
        confi = confint(f);%フィッティング係数の95%信頼区間(1行目下限, 2行目上限)
        passive_Ti(i) = 1.69e8*mass*(2*coeff(3)*sqrt(log(2))/lambda0)^2-Ti_instru_CH(i);
        passive_Timax(i) = 1.69e8*mass*(2*confi(2,3)*sqrt(log(2))/lambda0)^2-Ti_instru_CH(i);
        passive_Timin(i) = 1.69e8*mass*(2*confi(1,3)*sqrt(log(2))/lambda0)^2-Ti_instru_CH(i);
        passive_Em(i) = resolution*sum(spectra(:,i));
        if plot_CH_spectra
            idx_subp1 = mod(i-1,n_subp1)+1;%サブプロット位置番号
            if i == 1
                %figureの設定
                fig_num=2;
                figure(fig_num);
                gcf.Position=[0 500 400 300];
                sgtitle('Line Integrated Spectra')
            end
            subplot(col_subp1,raw_subp1,idx_subp1)
            if i <= n_subp1
                if i == 1
                    fit_x = transpose(linspace(lambda0-hw_plot_px*resolution,lambda0+hw_plot_px*resolution,100));
                    fit_y = feval(f,fit_x);
                    pfit1 = plot(fit_x,fit_y,'r-',lambda(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),spectra(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),'b+');%フィティングあり
                    pfit = repmat(pfit1,[1,1,n_subp1]);
                else
                    fit_y = feval(f,fit_x);
                    pfit(:,:,idx_subp1) = plot(fit_x,fit_y,'r-',lambda(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),spectra(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),'b+');%フィティングあり
                end
                xline(lambda0)
                title(['CH',num2str(ch(i))])
                xlabel('Wavelength [nm]')
                ylabel('Strength [a.u.]')
                xlim([lambda0-hw_plot_px*resolution lambda0+hw_plot_px*resolution])
                ylim([0 inf])
                legend('off')
            else
                fit_y = feval(f,fit_x);
                pfit(1,:,idx_subp1).YData = fit_y;
                pfit(2,:,idx_subp1).XData = lambda(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i);
                pfit(2,:,idx_subp1).YData = spectra(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i);
                title(['CH',num2str(ch(i))])
                drawnow
            end
        end
    end
end
%%
%%----------線積分温度、線積分発光強度二次元分布を計算------------
if cal_LineInt
    Ti_LineInt = zeros(270,1);%二次元温度[eV]
    Em_LineInt = zeros(270,1);%二次元発光強度[a.u.]
    %死んだCHを補間
    %ch_for_z:各zが(欠けがなければ)何チャンネルからなるのが理想かを示す
    ch_for_z=[16,16,28,30,30,30,30,30,28,16,16];
    ind=1;
    for i = 1:numel(z)
        if i < numel(z)
            Ti_LineInt(ind:ind+ch_for_z(i)-1) = pchip(ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1,passive_Ti(separation(i)+1:separation(i+1)),linspace(1,ch_for_z(i),ch_for_z(i)));
            Em_LineInt(ind:ind+ch_for_z(i)-1) = pchip(ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1,passive_Em(separation(i)+1:separation(i+1)),linspace(1,ch_for_z(i),ch_for_z(i)));
            ind=ind+ch_for_z(i);
        else
            Ti_LineInt(ind:ind+ch_for_z(i)-1) = pchip(ch(separation(i)+1:end)-ch(separation(i)+1)+1,passive_Ti(separation(i)+1:end),linspace(1,ch_for_z(i),ch_for_z(i)));
            Em_LineInt(ind:ind+ch_for_z(i)-1) = pchip(ch(separation(i)+1:end)-ch(separation(i)+1)+1,passive_Em(separation(i)+1:end),linspace(1,ch_for_z(i),ch_for_z(i)));
        end
    end
    %補間結果をプロットtodo
    if plot_LineInt_interp
        %figureの設定
        fig_num=3;
        figure(fig_num);
        gcf.Position=fig_position;
        % figure('Position',[1000 1000 600 200])
        tiledlayout(1,2)
        nexttile;
        p1 = plot(linspace(1,numel(p),numel(p)),Ti_LineInt(:,1),'ro-');
        hold on
        p2 = plot(ch(separation(1)+1:separation(1+1))-ch(separation(1)+1)+1,passive_Ti(separation(1)+1:separation(1+1)),'bo-');
        xlabel('CH')
        ylabel('Ti [eV]')
        xlim([1 numel(p)])
        title('Interpoltaion of Ti')
        legend(p1,'Interpolated')
        nexttile;
        p3 = plot(linspace(1,numel(p),numel(p)),Em_LineInt(:,1),'ro-');
        hold on
        p4 = plot(ch(separation(1)+1:separation(1+1))-ch(separation(1)+1)+1,passive_Em(separation(1)+1:separation(1+1)),'bo-');
        xlabel('CH')
        ylabel('Emission [a.u.]')
        xlim([1 numel(p)])
        title('Interpoltaion of Em')
        legend(p3,'Interpolated')
        hold off
        for i = 1:numel(z)
            sgtitle([num2str(i),'/',num2str(numel(z))])
            p1.YData = Ti_LineInt(:,i);
            p3.YData = Em_LineInt(:,i);
            if i < numel(z)
                p2.XData = ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1;
                p2.YData = passive_Ti(separation(i)+1:separation(i+1));
                p4.XData = p2.XData;
                p4.YData = passive_Em(separation(i)+1:separation(i+1));
                drawnow
            else
                p2.XData = ch(separation(i)+1:end)-ch(separation(i)+1)+1;
                p2.YData = passive_Ti(separation(i)+1:end);
                p4.XData = p2.XData;
                p4.YData = passive_Em(separation(i)+1:end);
                drawnow
            end
        end
    end
    %偶数列を上下反転(CHの順序が偶数列で反転するため)
%     for i = 1:numel(z)
%         if mod(i,2) == 0
%             Ti_LineInt(:,i) = flipud(Ti_LineInt(:,i));
%             Em_LineInt(:,i) = flipud(Em_LineInt(:,i));
%         end
%     end
end
%%
if plot_LineInt_result
    % figure('Position',[0 0 700 350])
    %figureの設定
    fig_num=4;
    figure(fig_num);
    gcf.Position=fig_position;
    
    tiledlayout(1,2)
    ax1 = nexttile;
    % [~,h] = contourf(z,p,Ti_LineInt,100);
    % h.LineStyle = 'none';
    % daspect([1 1 1])
    [xq,yq] = meshgrid(z,linspace(0.08,0.350,50));
    vq = griddata(C(:,7)*1e-3,C(:,8)*1e-3,Ti_LineInt,xq,yq);
    contourf(xq,yq,vq,50,'LineStyle','none');

    colormap(ax1,turbo)
    title('Line Integrated Ion Temperature')
    xlabel('Z [m]')
    ylabel('P [m]')
    c1 = colorbar;
    c1.Label.String = 'Ion Temperature [eV]';
    
    ax2=nexttile;
    [xq,yq] = meshgrid(z,linspace(0.08,0.350,50));
    vq = griddata(C(:,7)*1e-3,C(:,8)*1e-3,Em_LineInt,xq,yq);
    contourf(xq,yq,vq,50,'LineStyle','none');

    colormap(ax1,turbo)
    title('Line Integrated Ion Emission')
    xlabel('Z [m]')
    ylabel('P [m]')
    c2 = colorbar;
    c2.Label.String = 'Ion Emission [a.u.]';
end
%%
%%----------アーベル変換線積分温度、発光二次元分布を計算------------
if cal_2D
    %三角グリッドの補間によりアーベル変換に用いる2次元(λ,P)スペクトルを用意
    r = transpose(linspace(min(C(:,8))*1e-3,edge,num_r));
    dr = r(2) - r(1);
    lambdaA = linspace(-hw_lambdaA*resolution,hw_lambdaA*resolution,2*hw_lambdaA+1) + lambda0;
    lambdaA = transpose(lambdaA);
    spectra_interp = zeros(numel(lambdaA),numel(r),numel(z));
    for i = 1:numel(z)
        if i < numel(z)
            tri_x = zeros(2*hw_lambdaA+1,separation(i+1)-separation(i)+1);
            if mod(i-1,2) == 0
                for k = 1:separation(i+1)-separation(i)
                    idx_ch = separation(i)+k;
                    tri_x(:,k) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
                end
                tri_x(:,end) = tri_x(:,end-1);
                tri_y = [p(separation(i)+1:separation(i+1));edge];
                tri_z = [spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:separation(i+1)),zeros(numel(lambdaA),1)];
            else
                for k = 1:separation(i+1)-separation(i)
                    idx_ch = separation(i)+k;
                    tri_x(:,k+1) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
                end
                tri_x(:,1) = tri_x(:,2);
                tri_y = [edge;flipud(p(separation(i)+1:separation(i+1)))];
                tri_z = [zeros(numel(lambdaA),1),spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:separation(i+1))];
            end
        else
            tri_x = zeros(2*hw_lambdaA+1,numel(ch)-separation(i)+1);
            if mod(i-1,2) == 0
                for k = 1:numel(ch)-separation(i)
                    idx_ch = separation(i)+k;
                    tri_x(:,k) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
                end
                tri_x(:,end) = tri_x(:,end-1);
                tri_y = [p(separation(i)+1:end);edge];
                tri_z = [spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:end),zeros(numel(lambdaA),1)];
            else
                for k = 1:numel(ch)-separation(i)
                    idx_ch = separation(i)+k;
                    tri_x(:,k+1) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
                end
                tri_x(:,1) = tri_x(:,2);
                tri_y = [edge;flipud(p(separation(i)+1:end))];
                tri_z = [zeros(numel(lambdaA),1),spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:end)];
            end
        end
        buf_x = zeros(numel(tri_x),1);%CHごとの波長(CH数*ピクセル数)
        buf_y = zeros(numel(tri_x),1);%P
        buf_z = zeros(numel(tri_x),1);%スペクトル
        for j=1:numel(tri_x)
            index = mod(j-1,numel(tri_y));
            buf_x(j) = tri_x(ceil(j/numel(tri_y)),index+1);
            buf_y(j) = tri_y(index+1);
            buf_z(j) = tri_z(ceil(j/numel(tri_y)),index+1);
        end
        [grid_x, grid_y] = meshgrid(lambdaA,r);
        F = scatteredInterpolant(buf_x,buf_y,buf_z);
        grid_z = F(grid_x, grid_y);
        spectra_interp(:,:,i) = transpose(grid_z);
        % if plot_2D_interp%todo
            if i == 1
                figure('Position',[1000 0 500 300])
                [~,h] = contourf(grid_x,grid_y,grid_z);
                h.FaceAlpha = 0.7;
                h.LineStyle = 'none';
                colorbar
                colormap('jet')
                hold on
                T = delaunay(buf_x,buf_y);
                trp = triplot(T,buf_x,buf_y,'w');
                xlabel('Wavelength [nm]')
                ylabel('P [m]')
                xlim([min(lambdaA) max(lambdaA)])
                ylim([min(r) max(r)])
            else
                h.ZData = grid_z;
                T = delaunay(buf_x,buf_y);
                delete(trp)
                trp = triplot(T,buf_x,buf_y,'w');
                drawnow
            end
            title(['Interpolated Spectra at Z = ',num2str(z(i)), '[m]'])
        % end
    end
    %---アーベル変換----
    local_spectra = zeros(numel(lambdaA),numel(r),numel(z));
    spectra_interp = movmean(spectra_interp,5,1);%二次元スペクトルのスムージング(波長方向)
    spectra_interp = movmean(spectra_interp,round((numel(r)-1)/16),2);%二次元スペクトルのスムージング(P方向)
    for k=1:numel(z)
        derivative = diff(spectra_interp(:,:,k),1,2)/dr;
        for l=1:numel(lambdaA)
            for i=1:numel(r)
                for j=i:numel(r)-1
                    local_spectra(l,i,k) = local_spectra(l,i,k) - 1/pi*derivative(l,j)...
                        *log(r(j+1)*(1+sqrt(1-(r(i)/r(j+1))^2))/(r(j)*(1+sqrt(1-(r(i)/r(j))^2))));%Balandin's Abel inversion
                end
            end
        end
    end
    w_mov_l = 5;
    local_spectra = movmean(local_spectra,w_mov_l,1);%局所スペクトルのスムージング(波長方向)
    local_spectra = movmean(local_spectra,round((numel(r)-1)/16),2);%局所スペクトルのスムージング(R方向)
    Em_local = zeros(numel(r),numel(z));
    Ti_local = zeros(numel(r),numel(z));
    Ti_local_max = zeros(numel(r),numel(z));
    Ti_local_min = zeros(numel(r),numel(z));
    Ti_instru_local = zeros(numel(z),1);
    checker = ones(numel(r),numel(z));
    RSQ = zeros(numel(r),numel(z));
    for i=1:numel(z)
        if i <numel(z)
            %todo:確認
            % Ti_instru_local(i) = sum(Ti_instru_CH(separation(i)+1:separation(i+1)+1)/(separation(i+1)-separation(i)));
            Ti_instru_local(i) = sum(Ti_instru_CH(separation(i)+1:separation(i+1)-1)/(separation(i+1)-separation(i)));
        else
            Ti_instru_local(i) = sum(Ti_instru_CH(separation(i)+1:end)/(numel(ch)-separation(i)));
        end
    end
    %----フィッティング----
    close all;
    col_subp2 = 5;%サブプロット行数
    raw_subp2 = 9;%サブプロット列数
    n_subp2 = col_subp2*raw_subp2;
    col_subp_f = 2;%サブプロット行数
    raw_subp_f = 2;%サブプロット列数
    n_subp_f = col_subp_f*raw_subp_f;
    idx_f_fit = 1;
    col_subp_t = 2;%サブプロット行数
    raw_subp_t = 2;%サブプロット列数
    n_subp_t = col_subp_t*raw_subp_t;
    idx_t_fit = 1;
    for j=1:numel(z)
        for i=1:numel(r)
            idx_fit = (j-1)*numel(r)+i;%スペクトル番号
            input = local_spectra(:,i,j);
            for l=1:numel(lambdaA)
                if input(l) < 0
                    input(l) = -input(l) - min(abs(movmean(input,20)));
                end
            end
            try
                [f,gof] = fit(lambdaA(hw_lambdaA+1-hw_fit:hw_lambdaA+1+hw_fit),input(hw_lambdaA+1-hw_fit:hw_lambdaA+1+hw_fit),'gauss1');
                RSQ(i,j) = gof.rsquare;%決定係数(フィッティング精度指標)
                coeff = coeffvalues(f);%フィッティング係数
                % coeff(3)/resolution/sqrt(2);%IDLのcoeff[2]に等しい(確認用)
                confi = confint(f);%フィッティング係数の95%信頼区間(1行目下限, 2行目上限)
                Ti_local(i,j) = 1.69e8*mass*(2*coeff(3)*sqrt(log(2))/lambda0)^2-Ti_instru_local(j);
                Ti_local_max(i,j) = 1.69e8*mass*(2*confi(2,3)*sqrt(log(2))/lambda0)^2-Ti_instru_local(j);
                Ti_local_min(i,j) = 1.69e8*mass*(2*confi(1,3)*sqrt(log(2))/lambda0)^2-Ti_instru_local(j);
                Em_local(i,j) = resolution*sum(input);
                %外れ値(checker=0)の条件
                checker(i,j) = checker(i,j) * (abs(coeff(2)-lambda0) < 0.1);
                checker(i,j) = checker(i,j) * (coeff(1) > 0);
                checker(i,j) = checker(i,j) * (Em_local(i,j) > 0);
                checker(i,j) = checker(i,j) * (Ti_local(i,j) > 0);
                checker(i,j) = checker(i,j) * (RSQ(i,j) > 0.95);
                switch plot_2D_spectra
                    case 'off'
                    case 'all'                        
                        %todo:local spectraにチャンネル番号と温度表示する
                        idx_subp2 = mod(idx_fit-1,n_subp2)+1;%サブプロット位置番号
                        % if mod(idx_fit-1,n_subp2) == 0
                        if idx_fit-1== 0
                            figure('Position',fig_position);
                            sgtitle('Local Spectra at (Z,R) [mm]')
                        end
                        subplot(col_subp2,raw_subp2,idx_subp2)
                        if idx_fit <= n_subp2
                            if idx_fit == 1
                                fit_x = lambdaA;
                                fit_y = feval(f,fit_x);
                                pfit1 = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                                pfit = repmat(pfit1,[1,1,n_subp2]);
                            else
                                fit_y = feval(f,fit_x);
                                pfit(:,:,idx_subp2) = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                            end
                            xline(lambda0)
                            title(sprintf("(%.1f,%.1f), T_i=%.1f",z(j)*1e3,r(i)*1e3,Ti_local(i,j)))
                            xlabel('Wavelength [nm]')
                            ylabel('Strength [a.u.]')
                            xlim([min(lambdaA) max(lambdaA)])
                            ylim([0 inf])
                            legend('off')
                        else
                            fit_y = feval(f,fit_x);
                            pfit(1,:,idx_subp2).YData = fit_y;
                            pfit(2,:,idx_subp2).YData = input;
                            title(sprintf("(%.1f,%.1f), T_i=%.1f",z(j)*1e3,r(i)*1e3,Ti_local(i,j)))
                            drawnow
                        end
                    case 'good'
                        if checker(i,j) == 1
                            idx_subp_t = mod(idx_t_fit-1,n_subp_t)+1;%サブプロット位置番号
                            if idx_t_fit == 1
                                figure('Position',[400 500 400 300]);
                                sgtitle('Good Spectra at (Z,R) [mm]')
                            end
                            subplot(col_subp_t,raw_subp_t,idx_subp_t)
                            if idx_t_fit <= n_subp_t
                                if idx_t_fit == 1
                                    fit_x = lambdaA;
                                    fit_y = feval(f,fit_x);
                                    pfit_t1 = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                                    pfit_t = repmat(pfit_t1,[1,1,n_subp_t]);
                                else
                                    fit_y = feval(f,fit_x);
                                    pfit_t(:,:,idx_subp_t) = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                                end
                                title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                                xline(lambda0)
                                xlabel('Wavelength [nm]')
                                ylabel('Strength [a.u.]')
                                xlim([min(lambdaA) max(lambdaA)])
                                ylim([0 inf])
                                legend('off')
                            else
                                fit_y = feval(f,fit_x);
                                pfit_t(1,:,idx_subp_t).YData = fit_y;
                                pfit_t(2,:,idx_subp_t).YData = input;
                                title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                                drawnow
                            end
                            idx_t_fit = idx_t_fit + 1;
                        end
                    case 'bad'
                        if checker(i,j) == 0
                            idx_subp_f = mod(idx_f_fit-1,n_subp_f)+1;%サブプロット位置番号
                            if idx_f_fit == 1
                                figure('Position',[400 500 400 300])
                                sgtitle('Bad Spectra at (Z,R) [mm]')
                            end
                            subplot(col_subp_f,raw_subp_f,idx_subp_f)
                            if idx_f_fit <= n_subp_f
                                if idx_f_fit == 1
                                    fit_x = lambdaA;
                                    fit_y = feval(f,fit_x);
                                    pfit_f1 = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                                    pfit_f = repmat(pfit_f1,[1,1,n_subp_f]);
                                else
                                    fit_y = feval(f,fit_x);
                                    pfit_f(:,:,idx_subp_f) = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                                end
                                title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                                xline(lambda0)
                                xlabel('Wavelength [nm]')
                                ylabel('Strength [a.u.]')
                                xlim([min(lambdaA) max(lambdaA)])
                                ylim([0 inf])
                                legend('off')
                            else
                                fit_y = feval(f,fit_x);
                                pfit_f(1,:,idx_subp_f).YData = fit_y;
                                pfit_f(2,:,idx_subp_f).YData = input;
                                title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                                drawnow
                            end
                            idx_f_fit = idx_f_fit + 1;
                        end
                    otherwise
                        warning('Input error in plot_2D_spectra.')%ICCD.lineの入力エラー
                        return;
                end
            catch ME
                checker(i,j) = 0;
                warning(sprintf("Fitting failed in (Z,R) = (%.1f,%.1f)",z(j)*1e3,r(i)*1e3));
            end
        end
    end
    %NaNを除去
    for i=1:numel(r)
        for j=1:numel(z)
            if isnan(Ti_local(i,j))
                Ti_local(i,j) = 0;
                Ti_local_max(i,j) = 0;
                Ti_local_min(i,j) = 0;
                fprintf('Ti_2D(%d,%d) is NaN.',i,j);
            end
        end
    end
    %外れ値(checker=0)を補間
    Ti_local_smooth = filloutliers(Ti_local,"linear");
    Ti_local_smooth = movmean(Ti_local_smooth,5,1);
    Ti_local_smooth = movmean(Ti_local_smooth,round((numel(r)-1)/16),2);
    Ti_local_max_smooth = filloutliers(Ti_local_max,"linear");
    Ti_local_max_smooth = movmean(Ti_local_max_smooth,5,1);
    Ti_local_max_smooth = movmean(Ti_local_max_smooth,round((numel(r)-1)/16),2);
    Ti_local_min_smooth = filloutliers(Ti_local_min,"linear");
    Ti_local_min_smooth = movmean(Ti_local_min_smooth,5,1);
    Ti_local_min_smooth = movmean(Ti_local_min_smooth,round((numel(r)-1)/16),2);
    for i=1:numel(r)
        for j=1:numel(z)
            if checker(i,j) == 0
                Ti_local(i,j) = Ti_local_smooth(i,j);
                Ti_local_max(i,j) = Ti_local_max_smooth(i,j);
                Ti_local_min(i,j) = Ti_local_min_smooth(i,j);
                fprintf('Ti(%d,%d) was interpolated.\n',i,j);
            end
        end
    end
    Ti_local_offset = mean(Ti_local_min(end,:));
    Ti_local = Ti_local - Ti_local_offset;
    negative = find(Ti_local<0);
    Ti_local(negative) = zeros(size(negative));
end
%%
%%----------局所イオン温度のプロット
if plot_2D_result
    figure('Position',fig_position)
    tiledlayout(1,2)
    ax1 = nexttile;
    [~,h] = contourf(z,r,Ti_local,60);
    daspect([1 1 1])
    h.LineStyle = 'none';
    colormap(ax1,jet)
    clim([0,60]);
    title('Local Ion Temperature')
    xlabel('Z [m]')
    ylabel('R [m]')
    c1 = colorbar;
    c1.Label.String = 'Ion Temperature [eV]';
    ax2 = nexttile;
    [~,h] = contourf(z,r,Em_local,100);
    daspect([1 1 1])
    h.LineStyle = 'none';
    colormap(ax2,pink)
    title('Local Ion Emission')
    xlabel('Z [m]')
    ylabel('R [m]')
    c2 = colorbar;
    c2.Label.String = 'Ion Emission [a.u.]';
end
if plot_profile
    figure('Position',[1000 1000 600 450])
    colororder(hsv(numel(z)))
    for j=1:numel(z)
        leg = ['Z = ',num2str(z(j)),' [m]'];
        errorbar(Ti_local(:,j),r,(Ti_local_max(:,j)-Ti_local_min(:,j))/2,'horizontal','+-','LineWidth',1,'DisplayName',leg)
        hold on
    end
    legend
    title('Radial plofile of T_i in Z plane')
    xlabel('Ion Temperature [eV]')
    ylabel('R [m]')
    xlim([0 inf])
    ylim([min(r) max(r)])
end
%%
%%-----磁気面重ね書き
%close all
% clearvars -except date IDXlist doCheck
% addpath '/Users/rsomeya/Documents/GitHub/test-open'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス

%%%%%%%%%%%%%%%%%%%%%%%%
%280ch用新規pcbプローブのみでの磁気面（Bz）
%%%%%%%%%%%%%%%%%%%%%%%%

% run define_path.m
% %%%%%ここが各PCのパス
% %【※コードを使用する前に】環境変数定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
% pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
% pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
% pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
% pathname.save=getenv('savedata_path');%outputデータ保存先
% pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
% pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
% pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
% pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
% date=230714;
T=searchlog(T,node,date);
% IDXlist= 1; %[5:50 52:55 58:59];%[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];
n_data=numel(IDXlist);%計測データ数
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);

trange=400:800;%【input】計算時間範囲
n=40; %【input】rz方向のメッシュ数

%【input】プローブチェックか磁気面か
%doCheck = false;
%doCheck = true;

% figure('Position', [0 0 1500 1500],'visible','on');
for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    if shot == tfshot
        tfshot = [0,0];
    end
    i_EF=EFlist(i);
    TF=TFlist(i);
    if doCheck
        check_signal(date, shot, tfshot, pathname, n);
    else
        plot_psi200ch(date, shot, tfshot, pathname,n,i_EF,trange,IDXlist,z,r,Ti_local);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function plot_psi200ch(date, shot, tfshot, pathname, n,i_EF,trange,IDXlist,z,r,Ti_local)%230817加筆
% filename=strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% if exist(filename,"file")==0
%     disp('No rawdata file -- Start generating!')
%     rawdataPath = pathname.rawdata;
%     save_dtacq_data(dtacq_num, shot, tfshot,rawdataPath)
%     % return
% end
% load(filename,'rawdata');%1000×192
% 
% %正しくデータ取得できていない場合はreturn
% if numel(rawdata)< 500
%     return
% end

filename = strcat(pathname.processeddata,'/a039_',num2str(shot(1)),'.mat');
if exist(filename,'file') == 0
    doCalculation = true;
else
    doCalculation = false;
end


if doCalculation
%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));
C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
r_shift = 0.00;
ok = logical(C(:,14));
dtacq_num_list = C(:,1);
dtaq_ch = C(:,2);
polarity=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10)+r_shift;
ch=C(:,7);

if ismember(39,dtacq_num_list)
    filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    if exist(filename1,"file")==0
        disp('No rawdata file of a039 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
        % disp(['File:',filename1,' does not exit']);
        % return
    end
    a039_raw = importdata(filename1);
end
if ismember(40,dtacq_num_list)
    filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
    if exist(filename2,"file")==0
        disp('No rawdata file of a040 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
        % disp(['File:',filename2,' does not exit']);
        % return
    end
    a040_raw = importdata(filename2);
end

raw = zeros(1000,length(dtaq_ch));
for i = 1:length(dtaq_ch)
    if dtacq_num_list(i) == 39
        raw(:,i) = a039_raw(:,dtaq_ch(i));
    elseif dtacq_num_list(i) == 40
        raw(:,i) = a040_raw(:,dtaq_ch(i));
    end
end

b=raw.*coeff';%較正係数RC/NS
b=b.*polarity';%極性揃え

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=false(100,1);
ok_bt=ok_bz;
zpos_bz=zeros(100,1);
rpos_bz=zpos_bz;
zpos_bt=zpos_bz;
rpos_bt=zpos_bz;

%digital filter
windowSize = 8;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

for i=1:length(ch)
    b(:,i) = filter(bb,aa,b(:,i));
    b(:,i) = b(:,i) - mean(b(1:40,i));
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
        zpos_bz(ceil(ch(i)/2))=zpos(i);
        rpos_bz(ceil(ch(i)/2))=rpos(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
        zpos_bt(ceil(ch(i)/2))=zpos(i);
        rpos_bt(ceil(ch(i)/2))=rpos(i);
    end
end

% zprobepcb    = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
[zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
% [zq,rq]      = meshgrid(zprobepcb,rprobepcb);
[zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);
ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
for i = 1:length(ok_bt)
    if rpos_bt(i) > (r_shift)
        index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
        ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
    end
    index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
    ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
end

grid2D=struct(...
    'zq',zq,...
    'rq',rq,...
    'zprobepcb',zprobepcb,...
    'rprobepcb',rprobepcb,...
    'rprobepcb_t',rprobepcb_t,...1
    'ok_bz_matrix',ok_bz_matrix,...
    'ok_bt_matrix',ok_bt_matrix);
grid2D_probe = struct('zq',zq_probepcb,'rq',rq_probepcb,'rq_t',rprobepcb_t);

clear zq rq zprobepcb rprobepcb zq_probepcb rq_probepcb rprobepcb_t ok_bz_matrix ok_bt_matrix

% probecheck_script;

%data2Dcalc.m
r_EF   = 0.5 ;
n_EF   = 234. ;

if date<221119
    z1_EF   = 0.68;
    z2_EF   = -0.68;
else
    z1_EF   = 0.78;
    z2_EF   = -0.78;
end
[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
clear EF r_EF n_EF i_EF z_EF

data2D=struct(...
    'psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jr',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'trange',trange);

% ******************* no angle correction ********************
for i=1:size(trange,2)
    t=trange(i);

    %Bzの二次元補間(線形fit)
    vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
    B_z = -Bz_EF+vq;
    B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);

    % PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
    % data2D.psi(:,:,i) = flip(get_psi(flip(B_z,1),flip(grid2D.rq(:,1)),1),1);
    % このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bt(:,:,i)=B_t;
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end

else
    load(filename,'data2D','grid2D');
end
% ***********************************************

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

% プロット部分
%figure('Name',num2str(IDXlist),'Position', [0 0 1500 1500],'visible','on');
figure('Name',num2str(IDXlist),'Position', get(0,'screensize'),'visible','on');
start=50;%【変える】
dt = 2;

%  for m=1:50 %図示する時間
%      i=start+m.*dt; %end
%      t=trange(i);
%      if m == 1
%          [~,h1] = contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none');
%          hold on     
%          [~,h2] = contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1);
%          colormap(jet)
%          axis image
%          axis tight manual
%      else
%          h1.ZData = data2D.psi(:,:,i);
%          h2.ZData = squeeze(data2D.psi(:,:,i));
%          drawnow
%      end
%          title(string(t)+' us')
%          xlim([-0.17 0.17])
%      %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
% 
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
% %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
% %     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
% %     caxis([-0.01,0.01])%Bz
%      % clim([-0.1,0.1])%Bt
%     % clim([-5e-3,5e-3])%psi
% %     caxis([-500,400])%Et
% %     colorbar('Location','eastoutside')
%     %カラーバーのラベル付け
% %     c = colorbar;
% %     c.Label.String = 'Jt [A/m^{2}]';
%     % hold on
% %     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
% %     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
% %     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1)
% %     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
% %     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
%      % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
% %     xlabel('z [m]')
% %     ylabel('r [m]')
%  end

 %t_start=470+start;
 tate=4;%【変える】
 yoko=4;%【変える】
 maisu=tate*yoko;
 sgt=sgtitle(append("Date: ",num2str(date),", Shot: ",num2str(IDXlist)));
 sgt.FontSize=20;
 for m=1:maisu %図示する時間
     i=start+m.*dt; %end
     t=trange(i);
     subplot(tate,yoko,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
%     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
%     caxis([-0.01,0.01])%Bz
     % clim([-0.1,0.1])%Bt
    % clim([-5e-3,5e-3])%psi
%     caxis([-500,400])%Et
%     colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black','LineWidth',1)
%    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.1e-3:40e-3],'black','LineWidth',1)%[:ここがプロット間隔:]【たまに変える】
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
     % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
    hold off
    t=title(string(t)+' us');
    t.FontSize=15;
    %drawnow
    xlim([-0.17 0.17])
    ylim([0.1 0.3])
%     xlabel('z [m]')
%     ylabel('r [m]')
 end

%%----------イオン温度との重ね書き----------
target_t=468;%todo:自動取得
target_i=find(trange==target_t);
figure('Name',num2str(IDXlist),'Position', get(0,'screensize'),'visible','on');

%イオン温度
[~,h] = contourf(z,r,Ti_local,60);
daspect([1 1 1])
h.LineStyle = 'none';
colormap(jet)
clim([0,60]);
%title('Local Ion Temperature')
xlabel('Z [m]')
ylabel('R [m]')
c1 = colorbar;
c1.Label.String = 'Ion Temperature [eV]';
hold on;

contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,target_i)),50,'black','LineWidth',1);
hold off;
t=title(string(target_t)+' us');
t.FontSize=15;
%%----------

if doCalculation
    clearvars -except data2D grid2D shot pathname;
    filename = strcat(pathname.processeddata,'/a039_',num2str(shot(1)),'.mat');
    save(filename)
end

end

function save_dtacq_data(dtacq_num,shot,tfshot,rawdataPath)

% rawdataPath=getenv('rawdata_path'); %保存先

% %%%%実験オペレーションの取得
% DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
% T=getTS6log(DOCID);
% node='date';
% % pat=230707;
% pat=date;
% T=searchlog(T,node,pat);
% IDXlist=6;%[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];
% IDXlist = shot;
% n_data=numel(IDXlist);%計測データ数
% shotlist=T.a039(IDXlist);
% tfshotlist=T.a039_TF(IDXlist);
% EFlist=T.EF_A_(IDXlist);
% TFlist=T.TF_kV_(IDXlist);
% dtacqlist=39.*ones(n_data,1);

% dtacqlist=39;
% shotlist= 1819;%【input】dtacqの保存番号
% tfshotlist = 1817;

% dtacqlist=40;
% shotlist= 299;%【input】dtacqの保存番号
% tfshotlist = 297;

% date = 230706;%【input】計測日
% n=numel(shotlist);%計測データ数

%RC係数読み込み

% for i=1:n
%     dtacq_num=dtacqlist(i);
%     shot=shotlist(i);
%     tfshot=tfshotlist(i);
%     [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
%     save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
%     if tfshot>0
%         [rawdata]=getMDSdata(dtacq_num,shot,0);
%         % save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
%         save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
%     end
% end

% dtacq_num=dtacqlist(i);
% shot=shotlist(i);
% tfshot=tfshotlist(i);

[rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
if tfshot==0
    [rawdata]=getMDSdata(dtacq_num,shot,0);
    % save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
    save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
end

end

function check_signal(date, shot, tfshot, pathname,n)

sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));
C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
r_shift = 0.00;
ok = logical(C(:,14));
dtacq_num_list = C(:,1);
dtaq_ch = C(:,2);
polarity=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10)+r_shift;
ch=C(:,7);

if ismember(39,dtacq_num_list)
    filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    if exist(filename1,"file")==0
        disp('No rawdata file of a039 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
        % disp(['File:',filename1,' does not exit']);
        % return
    end
    a039_raw = importdata(filename1);
end
if ismember(40,dtacq_num_list)
    filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
    if exist(filename2,"file")==0
        disp('No rawdata file of a040 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
        % disp(['File:',filename2,' does not exit']);
        % return
    end
    a040_raw = importdata(filename2);
end

raw = zeros(1000,length(dtaq_ch));
for i = 1:length(dtaq_ch)
    if dtacq_num_list(i) == 39
        raw(:,i) = a039_raw(:,dtaq_ch(i));
    elseif dtacq_num_list(i) == 40
        raw(:,i) = a040_raw(:,dtaq_ch(i));
    end
end

b=raw.*coeff';%較正係数RC/NS
b=b.*polarity';%極性揃え
% b = smoothdata(b,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=false(100,1);
ok_bt=ok_bz;
zpos_bz=zeros(100,1);
rpos_bz=zpos_bz;
zpos_bt=zpos_bz;
rpos_bt=zpos_bz;

%digital filter
windowSize = 8;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

for i=1:length(ch)
    b(:,i) = filter(bb,aa,b(:,i));
    b(:,i) = b(:,i) - mean(b(1:40,i));
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
        zpos_bz(ceil(ch(i)/2))=zpos(i);
        rpos_bz(ceil(ch(i)/2))=rpos(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
        zpos_bt(ceil(ch(i)/2))=zpos(i);
        rpos_bt(ceil(ch(i)/2))=rpos(i);
    end
end

zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
[zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
for i = 1:length(ok_bt)
    if rpos_bt(i) > (r_shift)
        index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
        ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
    end
    index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
    ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
end

grid2D=struct(...
    'zq',zq,...
    'rq',rq,...
    'zprobepcb',zprobepcb,...
    'rprobepcb',rprobepcb,...
    'rprobepcb_t',rprobepcb_t,...1
    'ok_bz_matrix',ok_bz_matrix,...
    'ok_bt_matrix',ok_bt_matrix);

figure_switch = ["on","on","on","on","off","off"];%bz1,bz2,bt1,bt2,bz_vs_z,bt_vs_z

r = 7;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
y_upper_lim = 0.05;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.05;%縦軸プロット領域（b_z下限）
t_start=1;%横軸プロット領域（開始時間）
t_end=1000;%横軸プロット領域（終了時間）
t = 470;
% z_probe_pcb = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
z_probe_pcb = [-0.2975 -0.255 -0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17 0.255 0.2975];
n_z = length(z_probe_pcb);
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f1=figure(Visible=figure_switch(1));
f1.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j));
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe1-5') 

f2=figure(Visible=figure_switch(2));
f2.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe6-10')

f3=figure(Visible=figure_switch(3));
f3.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe1-5')

f4=figure(Visible=figure_switch(4));
f4.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe6-10')

% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close

% 横軸z, 縦軸Bzのプロット
f5=figure(Visible=figure_switch(5));
f5.WindowState = 'maximized';
styles = ["-*","-^","-v","-<","->","-o","-square","-diamond","-pentagram","-hexagram"];
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bz_zline=bz(t,zline);
    bz_zline(ok_bz(zline)==false)=NaN;
    plot(z_probe_pcb,bz_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('Before interpolation')

Bz_interped = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bz_zline=Bz_interped(zline);
    plot(linspace(min(zpos_bz),max(zpos_bz),n_z),bz_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('After interpolation')
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

% 横軸z, 縦軸Btのプロット
f6=figure(Visible=figure_switch(6));
f6.WindowState = 'maximized';
styles = ["-*","-^","-v","-<","->","-o","-square","-diamond","-pentagram","-hexagram"];
tiles = tiledlayout(2,1);
sgtitle(strcat('t=',num2str(t),' us'))
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bt_zline=bt(t,zline);
    bt_zline(ok_bt(zline)==false)=NaN;
    plot(z_probe_pcb,bt_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('Before interpolation')

Bt_interped = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
nexttile
hold on
for i=1:10
    zline=(1:10:n_z*10-9)+(i-1);
    bt_zline=Bt_interped(zline);
    plot(linspace(min(zpos_bz),max(zpos_bz),n_z),bt_zline,styles(i),'Color','k','MarkerSize',12)
    clear bz_zline
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
xline(z_probe_pcb,':','LineWidth',1.5)
legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')
title('After interpolation')
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

hidden = find(figure_switch == 'off');
figures = [f1,f2,f3,f4,f5,f6];
for i = hidden
    close(figures(i));
end

end
