%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 20/04/2022   %
%   Lab - experiment 2      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mass to time summary:
% Take damped harmonic data and fit it as:
% a*exp(b*x)*(cos(2*pi*x/c + 2*pi/d))+e
% Prints goodness of fit and plot

%% Parameters:
distance_error_range = 0.0172;
results_addr = 'D:\labs\simulation\lab_2\csv_files\5 G part 2.csv';

%% Grab lab results
results = readtable(results_addr);
y = results{:,2};
x = results{:,1};

%% Remove unwanted data

% Remove data captured before the experiment begins:
init_wave = find(y == min(y));
x = x-x(init_wave);
for i = 1:init_wave-1
	x(1) = [];
	y(1) = [];
end

%Remove data captured after the experiment ends:
act_zero = abs(y(length(y)))+distance_error_range;
for i = 60:length(x)
    if (abs(y(i))<act_zero) && (abs(y(i-60))<act_zero) 
        x(i:length(x)) = [];
        y(i:length(y)) = [];
        break
    end
end

%% Get fit parametes:

y = detrend(y);                                                                      % Remove Linear Trend
yu = max(y);
yl = min(y);
yr = (yu-yl);                                                                        % Range of  y’
yz = y-yu+(yr/2);
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                                 % Returns Approximate Zero-Crossing Indices Of Argument Vector
zt = x(zci(y));
per = 2*mean(diff(zt));                                                              % Estimate period
ym = mean(y);                                                                        % Estimate offset

init_fit = @(b,x)  b(1) .* exp(b(2).*x) .* (cos(2*pi*x./b(3) + 2*pi/b(4))) + b(5);   % Objective Function to fit
fcn = @(b) norm(init_fit(b,x) - y);                                                  % Least-Squares cost function
[s,nmrs] = fminsearch(fcn, [yr; -10;  per;  -1;  ym]);                               % Minimise Least-Squares

%% Fit:
fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', s');                % Use the parameters gathered as starting points.
f = fittype('a.*exp(b*x).*(cos(2*pi.*x/c + 2.*pi./d))+e','coefficients', {'a', 'b', 'c', 'd', 'e'}, 'options', fo);
final_fit = fit(x,y,f);

% Get goodness of fit:
fit_params = s'
fit_error_range = confint(final_fit)    

%% Plot:
figure
plot(x,y, '.')
hold on
plot(final_fit, '--r')
hold off
grid
box on
xlabel('Time(S)')
ylabel('Distance(M)')
legend('Original Data',  'Fitted Curve')
