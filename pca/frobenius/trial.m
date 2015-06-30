signal = zeros(1009, 91);

for i = 1:91
    signal(:,i) = data([ROIs.ndx{idx}], i);
end

ratio = zeros(1009, 90);
for i = 1:90
    ratio(:,i) = signal(:,i+1)./signal(:,i);
    ratio(isnan(ratio(:,i)),i) = 1;
end

fid = fopen('signal.txt','w');
for i = 1:1009
    for j = 1:91
        fprintf(fid, '%f ', signal(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

fid = fopen('ratio.txt','w');
for i = 1:1009
    for j = 1:90
        fprintf(fid, '%f ', ratio(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);