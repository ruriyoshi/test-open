function f = plot_multitime_spectrum(date,daybyday,BG,shot,ch,time,locate)
% ê‡ñæèëÇ≠ÇÃÇﬂÇÒÇ«Ç¢Ç©ÇÁÉTÉ{ÇËÇ‹Ç∑ÅB
    d = 1e6/133.6;
    alpha = 88;
    R = 998.8/25.4;
    h = 6.63*10^-34;
    c = 3.00*10^8;
    
    lambda = d*(sind(alpha) - sind(acosd((2.6+locate/25.4)/R)));
    
    
    data = [];
    data_BG = [];
    c = 1;
    cc = 1;
    for ii = daybyday
        for i = 1:ii
            data1 = get_one_ch(date(cc),BG(c),ch);
            data_BG = [data_BG data1];
            data1 = get_one_ch(date(cc),shot(c),ch);
            data = [data data1];
            c = c+1;
        end
        cc = cc+1;
    end
    
    
    f = figure;
    hold on
    for i = time
        txt = ['time=',num2str(i/10)];
        plot(lambda,data(i,:)-data_BG(i,:),'DisplayName',txt);
    end
    legend show
    grid on
end

function [Yhp, yhp] = lowpass_fft(data, pass)
    Y = fft(data);
    Fs = 1e6;       % Sampling frequency: 1Mhz
    T = 1/Fs;       % Sampling period
    L = 10000;      % Length of signal
    %t = (0:L-1)*T;  % Time vector
    NFFT = length(data);
    %Y = fft(data1);
    F = ((0:1/NFFT:1-1/NFFT)*Fs).';
    
    %{
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    figure;
    plot(f,P1) 
    %}
    
    Yhp = Y;
    %Yhp(F<=pass) = 0;
    Yhp(F>=pass & F<=Fs-pass) = 0;
    
    %{
    P2 = abs(Yhp/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    figure;
    plot(f,P1)
    %}
    
    yhp = ifft(Yhp,'symmetric');
    %plot(yhp);
end
