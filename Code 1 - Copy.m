%A hybrid inversion algorithm to obtain the resistivity of the uninvaded zone based on the array induction log using the folwing formation model Bayessian part
        % y = a+2(b-a)/1+exp(c*x)^d
        % c value is equal to 1 
clc; close all; clear all;

% Load the data(load your data on text filt )
load DD.txt;%replace the "DD with your uploaded dat"
data = DD; 
xdata1 = [0.254, 0.508, 0.762, 1.524, 2.286, 3.046];
resistivity = data(:, 3:8);% specify the  resistivity data column
shale = data(:, 2);%specify the shale volume data column it will be used for constarin
depth = data(:, 1);%specify the depth column
cutoff = 40; % cutoff value of shale volume
rr = resistivity;

[m, n] = size(rr);
ydata = zeros(m, n);
initial = zeros(m, 3);
xdata = xdata1;
% c = 1;
d = 12;

% Define the model function
model = @(x, xdata) x(1) + (2 * (x(2) - x(1))) ./ (1 + exp(xdata.^x(3)));

% Prepare initial guesses
for i = 1:m
    for j = 1:n
        ydata(i, j) = 1 ./ rr(i, j); % conductivity mS/m
    end
    initial(i, :) = [ydata(i, 6), ydata(i, 1), d];
end

num_guesses = size(initial, 1);

% Store posterior means for each initial guess
posterior_means = zeros(num_guesses, 3);

% Define likelihood function
likelihood = @(params, xdata, ydata) prod(normpdf(ydata, model(params, xdata), 0.1));

% Define prior distributions for parameters
prior_std = [0.1, 0.1,0.1]; % Prior standard deviation

% Initialize random number generator
rng('default');

for w = 1:num_guesses
    x0 = initial(w, :);
    niter = 10000; % Number of iterations
    chain = zeros(niter, 3); % Initialize parameter chain
    chain(1, :) = x0; % Initial parameter values
    accept_count = 0; % Counter for accepted samples

    % Perform sampling
    for k = 2:niter
        % Sample proposal from normal distribution
        proposal = chain(k - 1, :) + 0.01 * randn(1, 3) .* prior_std;
        % Compute acceptance probability
        alpha = min(1, likelihood(proposal, xdata, ydata(w, :)) / likelihood(chain(k - 1, :), xdata, ydata(w, :)));
        % Accept or reject proposal
        if rand < alpha
            chain(k, :) = proposal;
            accept_count = accept_count + 1;
        else
            chain(k, :) = chain(k - 1, :);
        end
    end

    % Discard burn-in samples
    burnin = 1000;
    chain = chain(burnin + 1:end, :);

    % Compute posterior means
    posterior_means(w, :) = mean(chain);
end

disp('Posterior means:');
disp(posterior_means);
 you can save your Posterior means data result on your computer as excell file

--------------------------------------------------------------------------------------------------------------------


