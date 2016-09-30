% root counter

function N = rootCounter(results)
N = sum(floor(abs(diff(sign(results)))/2),'omitnan');
end