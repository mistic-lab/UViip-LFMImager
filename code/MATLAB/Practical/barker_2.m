fs = 1e6
T = 494e-6
t = 0:1/fs:T; 

y = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];

barker = zeros(1, 494);

for m = 0:12
    for n = 1:38
        barker(m*38+n) = y(m+1);
    end
end

plot(barker)
        
audiowrite('./barker.wav', barker, fs);