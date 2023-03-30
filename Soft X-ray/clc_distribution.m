function EE = clc_distribution(M,K,U,s,v,VectorImage,plot_flag)

% M = N_projection;
% K = N_grid^2;
Z=U'*VectorImage.';

gmin=-15;gmax=15;
lg_gamma=gmin:1:gmax;
l_g = numel(lg_gamma);
gamma=10.^(lg_gamma);

% if K>M
%     v = v(:,1:M);
% end
% sigma = (diag(S)).';
% if M>K
%     sigma = [sigma zeros(1,M-K)];
% end

V1 = zeros(1,l_g);
V2 = zeros(1,l_g);
Vgamma = zeros(1,21);

for n=1:l_g
    rho = M*gamma(n)./(s.^2+M*gamma(n));
    V11 = rho.*(Z.');
    V1(n)=M*sum(V11.^2);
    V2(n)=(sum(rho))^2;
    Vgamma(n)=V1(n)/V2(n);
end

if plot_flag
    figure;plot(lg_gamma,Vgamma,'*');
    xlabel('logγ');
    ylabel('GCV');
end

[~,gamma_index]=min(Vgamma);
E = zeros(1,K);

for i=1:K
    if M>K
        v_1 = [v(i,:) zeros(1,M-K)];
    else
        v_1 = v(i,:);
    end
%     whos sigma
%     whos v_1
%     whos Z
    E1 = (s./(s.^2+M*10^(lg_gamma(gamma_index)))).*v_1.*(Z.');
    E(i)=sum(E1);
end

EE = reshape(E,sqrt(K),sqrt(K)); %ここで縦がr、横がzで左下が最小になる
% contourfで単調増加する軸から生成されたmeshgridを使ってプロットすると上下が反転する
EE = flipud(EE);

end