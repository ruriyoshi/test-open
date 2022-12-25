function [n_Amp,n_Ph] = toroidal_mode(t,x,signal)
t_start = t(1);
t_end = t(end);
len_freq = length(x);
signal = signal(t_start:t_end,:);
x_interp = 0:1/8:7/8;

% ********** main calculation **********
% https://i.stack.imgur.com/A2lb4.jpg
yshift_store = zeros(length(t),length(x));
for i = 1:length(t)
    sig = interp1(x,signal(i,:),x_interp,'makima','extrap' ); %補完関数
    y_fft = fft(sig);
    yshift_store(i,:) = fftshift(y_fft);
    %fshift = (-len_freq/2:len_freq/2)*((1/(abs(x(1)-x(2))))/len_freq);
end

% mode amplitude
n0_Amp = rot90(abs(yshift_store(:,len_freq/2+1)));
n1_Amp = 2*rot90(abs(yshift_store(:,len_freq/2)));
n2_Amp = 2*rot90(abs(yshift_store(:,len_freq/2-1)));
n3_Amp = 2*rot90(abs(yshift_store(:,len_freq/2-2)));

% mode phase
n1_Ph = -rot90(angle(yshift_store(:,len_freq/2)));
n2_Ph = -rot90(angle(yshift_store(:,len_freq/2-1)));
n3_Ph = -rot90(angle(yshift_store(:,len_freq/2-2)));

% return values
n_Amp = [n0_Amp;n1_Amp;n2_Amp;n3_Amp];
n_Ph = [n1_Ph;n2_Ph;n3_Ph];

end