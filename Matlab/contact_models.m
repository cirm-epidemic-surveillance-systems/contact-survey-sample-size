clear
close all

% Visualise some contact matrices constructed from the "continuous operator
% model" M(x,y) described in the Overleaf draft, with a continuous
% distribution of activity levels and different degrees of assortativity

% Model parameters

% Std. dev. in the log activity level distribution
Sigma = 0.1;

% Assortativity constant (in additon to Epsilon=0 and Epsulon=1 which will always
% be plotted)
Epsilon = 0.3;

% Width of assortativity kernel
ak = [0.1, 0.5];
nCases = length(ak);

% Number of bins to discretise and plot the contact matrix
nBins = 100;




% Create a grid of nBins x values placed at the midpoint of each "bin"
dx = 1/nBins;
x = dx/2:dx:(1-dx/2);
nx = length(x);

% Calculate v as the inverse CDF of the log normal activity level distribution
v = logninv(x, 0, Sigma);

% Calculate approximate moments of the distribution
Ev = dx*sum(v);

% Define matrices of x and y values for calculating M(x,y)
[X, Y] = meshgrid(x, x);


% Proportionate mixing model
M_PM = v'.*v/Ev;



% Initialise stack of matrices M_AM
M_AM = zeros(nx, nx, nCases);

% Calculate the fully assortative mixing model (epislon=1) for each value of ak
for iCase = 1:nCases
    % Calculate the kernel as a fuction of Y-X
    gk = exp(-(Y-X).^2/(2*ak(iCase)^2));

    % Calculate the produce v(y)*g(y-x)
    C = v'.*gk;

    % Calculte denominator of equation for M
    den = dx*sum(C, 1);

    % Assortative mixing matrix
    M_AM(:, :, iCase) = v.*C./den;
end


% Calculate the aggregate contacts for someone as a function of their
% activity class by summing columns of the matrix
aggCont_PM = dx*sum(M_PM, 1);
aggCont_AM = squeeze(dx*sum(M_AM, 1))';


figure(1)
tiledlayout(2, 3, "TileSpacing", "compact")
tileInd = [1, 2, 3, 5, 6];
for iMatrix = 1:5
    nexttile(tileInd(iMatrix));
    iCase = mod(iMatrix, 2)+1;
    if iMatrix == 1
        imagesc(x, x, M_PM);
        title('prop. mixing')
    elseif iMatrix == 2 | iMatrix == 3
        imagesc(x, x, (1-Epsilon)*M_PM + Epsilon*M_AM(:, :, iCase));
        title(sprintf('ass. mixing with eps = %.1f, a = %.1f', Epsilon, ak(iCase)))
    elseif iMatrix == 4 | iMatrix == 5
        imagesc(x, x, M_AM(:, :, iCase));
        title(sprintf('ass. mixing with eps = 1, a = %.1f', ak(iCase)))
    end
    h = gca;
    h.YDir = 'normal';
    colorbar;
    xlabel('activity level quantile of individual (x)')
    ylabel('activity level quantile of contact (y)')
end


figure(2)
plot(x, aggCont_PM)
hold on
plot(x, aggCont_AM)
plot(x, v, '--')
grid on
xlabel('activity level quantile')
ylabel('total contact rate')
legend('PM', 'AM1', 'AM2', 'log-norm dist', 'location', 'northwest')
title('check total contact rates match the specified distribution')
