function [vectorized,sim_index] = cosine_similarity(Spikes,bin)
x = length(Spikes(:,1));
y =length(Spikes(1,:));
try
    a = zeros(x,y/bin);
catch ME
    if rem(y,bin) ~= 0
        bin = bin+rem(y,bin);
        warning('Bin size does not match vector. Assigning new bin value');
    else
        msg = ('Bin size do not match vector length');
        causeException = MException('MATLAB:myCode:dimensions',msg);
        ME = addCause(ME,causeException);
        rethrow(ME);
    end
end
a = zeros(x,y/bin);
for i = 1:y/bin
    position = bin*(i-1)+1;
    for ii=1:bin
        a(position+(ii-1),i) = 1;
    end 
end
vectorized = Spikes*a;
timei = length(vectorized(1,:));
sim_index = zeros(timei,timei);
for j = 1:timei
    for n = 1:timei
        sim_index(j,n) = dot(vectorized(:,n),vectorized(:,j))/(norm(vectorized(:,n))*norm(vectorized(:,j)));
        if isnan(sim_index(j,n))
            sim_index(j,n) = 0;
        end
    end
end
% h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
end
