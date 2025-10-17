%Create a parallel cluster object
c = parcluster('local'); % local is a cluster profile name
%sz = str2num([getenv('SLURM_CPUS_PER_TASK')]); 
%if isempty(sz), sz = maxNumCompThreads; end
%if isempty(gcp('nocreate')), c.parpool(sz); end

p = c.parpool(10);

primeNumbers = primes(uint64(2^20));
compositeNumbers = primeNumbers.*primeNumbers(randperm(numel(primeNumbers)));
factors = zeros(numel(primeNumbers),2);

tic;
parfor idx = 1:numel(compositeNumbers)
    factors(idx,:) = factor(compositeNumbers(idx));
end
time=toc;

%Display the time of computation
fprintf('The parallel method is executed in seconds=')
disp(time)

delete(gcp('nocreate'))