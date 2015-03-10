frames2 = frames';
[n m] = size(frames2);
figure

for i = 1:m
    plot(frames(i,:))
    hold on
    keyboard
end