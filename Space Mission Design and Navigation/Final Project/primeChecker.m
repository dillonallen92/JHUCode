clear, clc, close all;


primeFunction(4)

function isPrime = primeFunction(testVal)
    isPrime = primeChecker(testVal, testVal-1);
end

function primeCheck = primeChecker(testVal, divVal)
    if (testVal == 0 || testVal == 1)
        primeCheck = 0;
    end
    
    if(mod(testVal, divVal) == 0)
        primeCheck = 0;
    end
    
    if (testVal == 2)
        primeCheck = 1;
    elseif (mod(testVal, 1) == 0 && mod(testVal,testVal) == 0 && mod(testVal, divVal) ~=0)
        primeCheck = 1;
    else
        primeChecker(testVal, divVal-1);
    end
        
end