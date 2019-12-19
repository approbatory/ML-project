% Demo for a simple progress bar without known total number of iterations
%
% Author:  J.-A. Adrian (JA) <jensalrik.adrian AT gmail.com>
% Date  :  21-Jun-2016 17:13:27
%


addpath('..');

numIterations = 50;



%% Simple setup WITHOUT known number of iterations

obj = ProgressBar();

for iIteration = 1:numIterations
    pause(0.1);
    
    obj.step(1, [], []);
end
obj.release();




%% Simple setup WITHOUT known number of iterations and with custom title

obj = ProgressBar([], 'Title', 'Test');

for iIteration = 1:numIterations
    pause(0.1);
    
    obj.step(1, [], []);
end
obj.release();





% End of file: a_ProgressBarWithoutTotal_demo.m
