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