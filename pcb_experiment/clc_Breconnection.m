function B_recconection = clc_Breconnection(grid2D,data2D)
z_space = grid2D.zq(1,:);
r_space = grid2D.rq(:,1);
trange = data2D.trange;
n = numel(z_space);
B_recconection = zeros(numel(trange),1);

frame=numel(trange);
%各時間、各列(z)ごとのpsiの最大値
[max_psi,max_psi_r]=max(data2D.psi,[],1);
max_psi=squeeze(max_psi);

for i=1:frame
    r_ind=max_psi_r(:,:,i);
    max_r = zeros(numel(r_ind),1);
    for j = 1:numel(r_ind)
        max_r(j) = r_space(r_ind(j));
    end
%     max_psi_ind=find(islocalmax(smooth(max_psi(:,i)),'MaxNumExtrema', 2));
%     max_psi_ind_left=find(max(smooth(max_psi(1:n/2,i))));
%     max_psi_ind_right=find(max(smooth(max_psi(n/2:n,i))));

    if numel(find(islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1)))==0
        continue
    else
        min_psi_ind=islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1);
        xr=r_ind(min_psi_ind);
        if xr==1 || xr==n %r両端の場合は検知しない
            continue
        else
%             figure;plot(z_space,max_psi(:,i));
            pp = spline(z_space,max_psi(:,i));
            p_der = fnder(pp,1);
%             A = (-1./(2*pi*max_r));
%             whos A
%             B = ppval(p_der,z_space);
%             whos B
            B_r = (-1./(2*pi*max_r)).'.*ppval(p_der,z_space);
%             whos B_r
%             figure;plot(z_space,B_r);
%             B_rec
%             break
            B_recconection(i) = (abs(min(B_r))+max(B_r))/2;
        end
    end

end
end
