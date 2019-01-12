function printMatrix(M)
for k=1:length(M)
    fprintf('%10.4f ',M(k,:));
    fprintf('\n')
end
end