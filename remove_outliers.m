function res = remove_outliers(data, dev)
    mu = mean(data,1);
    sigma = std(data,1,1);
    tmp = abs(data-mu);
    logIdx = all([tmp(:,1) <= dev*sigma(1), tmp(:,2) <= dev*sigma(2)],2);
    indices = 1:length(logIdx);
    indices = indices(logIdx);

    res = data(indices,:);
end