function assignment1_random_generators()
%% CH5030 Molecular Thermodynamics - Assignment 1
% Random number generation and statistical tests.
%CONNECT ME ON WHATSAPP 9711603123
% The program constructs all numerical data internally. No external dataset
% is needed.

clc; close all;

nValues = 10000;
seedLCG = 13579;
seedBits = makeSeedBits(400, 2026);

streams(1) = makeStream( ...
    'LCG: a = 48271, b = 0, m = 2^31 - 1', ...
    lcgStream(nValues, seedLCG, 48271, 0, 2^31 - 1));

streams(2) = makeStream( ...
    'LCG: a = 1664525, b = 1013904223, m = 2^32', ...
    lcgStream(nValues, seedLCG, 1664525, 1013904223, 2^32));

streams(3) = makeStream( ...
    'Shift register: L = 16, p = 103, q = 250', ...
    shiftRegisterStream(nValues, seedBits, 16, 103, 250));

streams(4) = makeStream( ...
    'Shift register: L = 32, p = 103, q = 250', ...
    shiftRegisterStream(nValues, seedBits, 32, 103, 250));

rng(2026, 'twister');
streams(5) = makeStream('Mersenne Twister: MATLAB rand', rand(nValues, 1));

fprintf('\nCH5030 Assignment 1: Random number tests\n');
fprintf('Number of values in each stream: %d\n\n', nValues);
fprintf('%-48s %10s %12s %12s %12s\n', ...
    'Generator', 'Mean', 'Mean z', 'Runs z', 'Runs');
fprintf('%s\n', repmat('-', 1, 99));

for i = 1:numel(streams)
    u = streams(i).values;
    meanStats = meanTest(u);
    runStats = runsTestAboutMedian(u, 0.5);
    streams(i).meanStats = meanStats;
    streams(i).runStats = runStats;

    fprintf('%-48s %10.5f %12.3f %12.3f %12d\n', ...
        streams(i).name, meanStats.sampleMean, meanStats.zScore, ...
        runStats.zScore, runStats.numberOfRuns);
end

plotHistograms(streams);
plotEmpiricalCdfFits(streams);

fprintf('\nInterpretation guide:\n');
fprintf('  For the mean test, z-scores close to 0 indicate a sample mean near 0.5.\n');
fprintf('  For the runs test, z-scores close to 0 indicate no obvious long-run clustering around 0.5.\n');
fprintf('  The Mersenne Twister stream is included as the reference generator requested in the assignment.\n');
end

%% Local functions

function item = makeStream(name, values)
    item.name = name;
    item.values = values(:);
    item.meanStats = [];
    item.runStats = [];
end

function u = lcgStream(n, seed, a, b, m)
    states = zeros(n, 1);
    state = uint64(mod(seed, m));
    multiplier = uint64(a);
    increment = uint64(b);
    modulus = uint64(m);

    for k = 1:n
        state = mod(multiplier * state + increment, modulus);
        states(k) = double(state);
    end

    u = states ./ double(modulus);
end

function bits = makeSeedBits(nBits, seed)
    rng(seed, 'twister');
    bits = rand(nBits, 1) > 0.5;

    if all(bits == 0)
        bits(1) = true;
    end
end

function u = shiftRegisterStream(nValues, seedBits, wordLength, p, q)
    if numel(seedBits) <= q
        error('The seed bit sequence must be longer than q.');
    end

    seedLen = numel(seedBits);
    totalBits = nValues * wordLength;
    bits = false(seedLen + totalBits, 1);
    bits(1:seedLen) = logical(seedBits(:));

    for k = (seedLen + 1):(seedLen + totalBits)
        bits(k) = xor(bits(k - p), bits(k - q));
    end

    generated = bits((seedLen + 1):(seedLen + totalBits));
    words = reshape(generated, wordLength, nValues).';
    weights = 2 .^ ((wordLength - 1):-1:0);
    integers = double(words) * weights(:);

    u = integers ./ 2^wordLength;
end

function stats = meanTest(u)
    n = numel(u);
    expectedMean = 0.5;
    expectedVariance = 1 / 12;
    standardError = sqrt(expectedVariance / n);

    stats.sampleMean = mean(u);
    stats.zScore = (stats.sampleMean - expectedMean) / standardError;
end

function stats = runsTestAboutMedian(u, threshold)
    labels = u >= threshold;
    labels = labels(:);

    n1 = sum(labels);
    n2 = numel(labels) - n1;
    runs = 1 + sum(labels(2:end) ~= labels(1:end-1));

    expectedRuns = 1 + 2 * n1 * n2 / (n1 + n2);
    varianceRuns = (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / ...
        ((n1 + n2)^2 * (n1 + n2 - 1));

    stats.numberOfRuns = runs;
    stats.expectedRuns = expectedRuns;
    stats.zScore = (runs - expectedRuns) / sqrt(varianceRuns);
end

function plotHistograms(streams)
    fig = figure('Name', 'Assignment 1 - Histograms', 'Color', 'w');
    set(fig, 'Position', [100 100 1400 1000]);

    for i = 1:numel(streams)
        subplot(3, 2, i);
        [counts, centers] = hist(streams(i).values, 30);
        binWidth = centers(2) - centers(1);
        density = counts ./ (sum(counts) * binWidth);
        bar(centers, density, 1.0, 'FaceColor', [0.15 0.45 0.70], ...
            'EdgeColor', 'none');
        hold on;
        plot([0 1], [1 1], 'k--', 'LineWidth', 1);
        title(streams(i).name, 'Interpreter', 'none', 'FontSize', 9);
        xlabel('u');
        ylabel('Density');
        xlim([0 1]);
        grid on;
    end

    saveFigure(fig, 'assignment1_histograms.png');
end

function plotEmpiricalCdfFits(streams)
    fig = figure('Name', 'Assignment 1 - Empirical CDF fits', 'Color', 'w');
    set(fig, 'Position', [100 100 1400 1000]);

    for i = 1:numel(streams)
        subplot(3, 2, i);
        x = sort(streams(i).values);
        y = ((1:numel(x)).' - 0.5) ./ numel(x);

        plot(x, y, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 4);
        hold on;

        if exist('fit', 'file') == 2
            fittedLine = fit(x, y, 'poly1');
            yFit = fittedLine.p1 .* x + fittedLine.p2;
            plot(x, yFit, 'r-', 'LineWidth', 1.5);
            fitLabel = sprintf('Curve fit: y = %.3fx %+.3f', ...
                fittedLine.p1, fittedLine.p2);
        else
            coeff = polyfit(x, y, 1);
            yFit = polyval(coeff, x);
            plot(x, yFit, 'r-', 'LineWidth', 1.5);
            fitLabel = sprintf('polyfit fallback: y = %.3fx %+.3f', ...
                coeff(1), coeff(2));
        end

        plot([0 1], [0 1], 'k--', 'LineWidth', 1);
        title(streams(i).name, 'Interpreter', 'none', 'FontSize', 9);
        xlabel('Sorted random number');
        ylabel('Empirical CDF');
        legend('Generated data', fitLabel, 'Ideal uniform CDF', ...
            'Location', 'southeast', 'Interpreter', 'none');
        xlim([0 1]);
        ylim([0 1]);
        grid on;
    end

    saveFigure(fig, 'assignment1_empirical_cdf_fits.png');
end

function saveFigure(fig, fileName)
    try
        print(fig, fileName, '-dpng', '-r300');
    catch
        saveas(fig, fileName);
    end
end
