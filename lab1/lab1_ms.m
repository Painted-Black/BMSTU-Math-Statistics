function lab1_ms()
    %clear all;
    X = readFromFile('C:\Users\novoc\Desktop\6sem\data.csv');
    X = sort(X); % Âàðèàöèîííûé ðÿä
    Xmin = X(1);
    Xmax = X(end);
    R = Xmax - Xmin;
    n = length(X);
    
    mu = sum(X) / n;
    ssqr = getCorrectedSampleVariance(X, mu);
    
    p = countSubintervals(n);
    
    fprintf('Ìèíèìàëüíîå çíà÷åíèå âûáîðêè: %.6f\n', Xmin);
    fprintf('Ìàêñèìàëüíîå çíà÷åíèå âûáîðêè: %.6f\n', Xmax);
    fprintf('Ðàçìàõ âûáîðêè: %.6f\n', R);
    fprintf('Ðàçìåð âûáîðêè: %d\n', n);
    fprintf('Âûáîðî÷íîå ìàòåìàòè÷åñêîå îæèäàíèå: %.6f\n', mu);
    fprintf('Èñïðàâëåííàÿ âûáîðî÷íàÿ äèñïåðñèÿ: %.6f\n', ssqr);
    fprintf('Âñåãî èíòåðâàëîâ: %d\n', p);
    
    [J, count] = group(X, p);
    drawHist(X, J, count);
    hold on; % hold on ñîõðàíÿåò ãðàôèêè â òåêóùåé ñèñòåìå êîîðäèíàò òàê, 
             % ÷òîáû íîâûå ãðàôèêè, äîáàâëåííûå ê îñÿì, íå óäàëÿëè ñóùåñòâóþùèå ãðàôèêè
    f(X, mu, ssqr, p);
    
    figure; % figure ñîçäàåò íîâîå îêíî ðèñóíêà
    drawDist(X);
    hold on;
    F(X, mu, ssqr, p);    
end

function F(X, MX, DX, p)
    R = X(end) - X(1);
    delta = R / p;
    Xn = (MX - R) : delta / 20 : (MX + R);
    Y = 0.5 * (1 + erf((Xn - MX) / sqrt(2 * DX)));
    plot(Xn, Y, 'r');
end

function drawDist(sample)
    [f, x] = ecdf(sample); % âîçâðàùàåò ýìïèðè÷åñêóþ èíòåãðàëüíóþ ôóíêöèþ ðàñïðåäåëåíèÿ (cdf) f,
                           % âû÷èñëåííóþ â òî÷êàõ x, èñïîëüçóÿ äàííûå â âåêòîðå y.
    stairs(x, f), grid;
end

function f (X, MX, DX, m)
    R = X(end) - X(1);
    delta = R / m;
    sigma = sqrt(DX);
    Xn = (MX - R) : delta / 20 : (MX + R);
    Y = normpdf(Xn, MX, sigma);
    plot(Xn, Y);
end

function drawHist(sample, J, count)
    p = length(count);
    n = length(sample);
    delta = (sample(n) - sample(1)) / p;
    xes = zeros(1, p + 1);
    xes(1) = 0;
    for i = 2 : p
        xes(i) = count(i) / (n * delta);
    end
    stairs(J, xes), grid;
end

function [J, count] = group(sample, p)
    delta = (sample(end) - sample(1)) / p;
    count = zeros(1, p);
    
    J = sample(1):delta:sample(end);
    
    fprintf('Èíòåðâàëû:\n');
    for i = 1 : p - 1
        fprintf('[%.6f; %.6f), ', J(i), J(i+1));
    end
    fprintf('[%.6f; %.6f]\n ', J(end - 1), J(end));
    
    for i = 1 : length(sample)
        cur = sample(i);
        for j = 1 : p - 1
            if ((cur >= J(j)) && (cur < J(j + 1)))
                count(j) = count(j) + 1;
            end
        end
        if ((cur >= J(end - 1)) && (cur <= J(end)))
            count(end) = count(end) + 1;
        end
    end
end

function p = countSubintervals(sample_size)
    p = floor(log2(sample_size)) + 1;
end
    
function ssqr = getCorrectedSampleVariance(sample, sample_mean)
    ssqr = sum((sample - sample_mean).^2) / (length(sample) - 1);
end

function X = readFromFile(filename)
    X = readmatrix(filename);
end
