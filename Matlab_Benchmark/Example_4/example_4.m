primeNumbers = primes(uint64(2^20));
compositeNumbers = primeNumbers.*primeNumbers(randperm(numel(primeNumbers)));
factors = zeros(numel(primeNumbers),2);

tic;
for idx = 1:numel(compositeNumbers)
    factors(idx,:) = factor(compositeNumbers(idx));
end
time=toc;

%Display the time of computation
fprintf('The serial method is executed in seconds=')
disp(time)



