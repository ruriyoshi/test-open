function anime_fluctuation(data)

for i = 1:24
    figure;
    ylim([-2000 2000]);
    plot(data(i,:));
    pause(0.2);
    close
    if i == 24
        break
    end
end
end