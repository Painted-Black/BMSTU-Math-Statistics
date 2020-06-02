function lab2()
    X = [-2.79,-3.01,-4.07,-2.85,-2.43, -3.20,-3.72,-4.27,-5.48,-2.38, -4.69, ...
        -4.34,-5.08,-5.01,-4.08, -4.20,-4.74,-1.88,-3.25,-2.78, -3.56,-3.54, ...
        -3.79,-3.18,-5.08,-4.30,-2.86,-2.45,-3.08,-3.22,-2.76,-3.20,-3.33, ...
        -4.91,-4.06,-3.81,-3.96,-3.65,-3.77,-4.60,-5.21,-2.67,-1.95,-2.43, ...
        -1.73,-2.50,-3.96,-3.75,-2.70,-4.26,-3.42,-4.07,-4.74,-3.00,-4.37, ...
        -5.42,-5.00,-4.08,-2.46,-4.33,-4.08,-3.72,-4.09,-2.96,-3.71,-1.51, ...
        -3.70,-6.48,-4.26,-4.39,-3.16,-4.63,-2.66,-2.22,-4.79,-2.46,-3.69, ...
        -3.35,-2.32,-4.17,-3.85,-4.93,-2.05,-3.15,-3.49,-5.70,-2.53,-3.85, ...
        -4.32,-3.37,-3.98,-3.74,-5.28,-2.56,-3.21,-3.10,-3.78,-3.36,-3.32, ...
        -2.59,-2.45,-3.34,-3.20,-4.14,-4.00,-4.79,-4.02,-4.58,-4.45,-3.69, ...
        -4.53,-3.98,-4.51,-4.44,-3.78,-4.24,-4.00,-2.46,-2.58,-4.04];
    
    N = length(X);
    
    mu = sExpectation(X);
    fprintf('mu = %.6f\n', mu);
    
    s_2 = correctedSampleVariance(X);
    fprintf('s_2 = %.6f\n', s_2);
    
    gamma = 0.9;
    alpha = (1.0 - gamma) / 2.0;
    fprintf('gamma = %.2f, apha = %.6f, N = %d\n', gamma, alpha, N);
    
    [lmu, umu] = getMXBorders(gamma, s_2, mu, N);
    fprintf('Нижняя гамма-доверительная граница для мат. ож.(x_N) = %.6f\n', lmu);
    fprintf('Верхняя гамма-доверительная граница для мат. ож.(x_N) = %.6f\n', umu);
    
    [ls, hs] = getDXBorders(gamma, s_2, N);
    fprintf('Нижняя гамма-доверительная граница для дисперсии (x_N) = %.6f\n', ls);
    fprintf('Верхняя гамма-доверительная граница для дисперсии (x_N) = %.6f\n', hs);
    
    figure(1);
    grid on;
    hold on;
    xlabel('n');
    ylabel('\mu');
    graphMX(X, N, gamma);
    
    figure(2);
    grid on;
    hold on;
    xlabel('n');
    ylabel('\sigma');
    graphDX(X, N, gamma);
end

function graphDX(X, n, gamma)
    s2s = zeros(n, 1);
    lowerSigma = zeros(n, 1);
    upperSigma = zeros(n, 1);
    
    for i = 1:n
        currentSample = X(1:i);
        [s2s(i)] = correctedSampleVariance(currentSample);
        [lowerSigma(i), upperSigma(i)] = getDXBorders(gamma, s2s(i), i);
    end
    
    plot([1, n], [s2s(n), s2s(n)], 'g');
    plot(lowerSigma, 'b');
    plot(upperSigma, 'r');
    plot(s2s, 'k');
    legend('S^2(x_N)', '_{--}\sigma^2(x_n)', '^{--}\sigma^2(x_n)', 'S^2(x_n)');
end

function [ls, hs] = getDXBorders(gamma, s_2, n)
    % неизвестны матожидание и дисперсия, оцениваем дисперсию;
    % статистика ~chi2(n-1)

    alpha1 = (1 + gamma) / 2;
    alpha2 = (1 - gamma) / 2;
    
    quantile1 = chi2inv(alpha1, n - 1);
    quantile2 = chi2inv(alpha2, n - 1);
    
    ls = ((n - 1) * s_2) / quantile1;
    hs = ((n - 1) * s_2) / quantile2;
end

function graphMX(X, n, gamma)
    mus = zeros(n, 1);
    s2s = zeros(n, 1);
    lowerMus = zeros(n, 1);
    upperMus = zeros(n, 1);
    
    for i = 1:n
        currentSample = X(1:i);
        [mus(i)] = sExpectation(currentSample);
        [s2s(i)] = correctedSampleVariance(currentSample);
        [lowerMus(i), upperMus(i)] = getMXBorders(gamma, s2s(i), mus(i), i);
    end
    plot([1, n], [mus(n), mus(n)], 'g');
    plot(lowerMus, 'b');
    plot(upperMus, 'r');
    plot(mus, 'k');
    legend('\mu\^(x_N)', '_{--}\mu^(x_n)', '^{--}\mu^(x_n)', '\mu\^(x_n)');
end

function [lm, hm] = getMXBorders(gamma, s_2, mu, n)
    % неизвестны мат. ожидание и дисперсия, оцениваем матожидание;
    % статистика ~St(n-1)
    alpha = (1.0 + gamma) / 2.0; % alpha1 = alpha2
    
    quantile = tinv(alpha, n - 1); % расчет значений квантили распр-я Стьюдента 
                                   % для значений вероятности alpha и степени свободы n - 1.
    border = (sqrt(s_2) * quantile) / sqrt(n);
    
    lm = mu - border;
    hm = mu + border;
end

function s_2 = correctedSampleVariance(X)
    s_2 = var(X); % исправленная выборочная дисперсия
end

function mu = sExpectation(X)
    mu = mean(X); % mean возвращает арифметическое среднее значение элементов массива
                  % выборочное мат. ожидание
end