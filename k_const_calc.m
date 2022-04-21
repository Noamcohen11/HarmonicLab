%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Noam Cohen 29/03/2022   %
%   Lab - experiment 2      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% K const test:

%%%%%%%%% Constants %%%%%%%%%

Weights = [
0
4.8
9.8
19.9
50.2
]./1000;

Weights_error = [
0
0.05
0.05
0.1
0.1
]./1000;


X = [0.172
0.086
0.275
0.361
0.447
0.636
0.911
0.997
1.186
1.1
1.548
1.272
1.376
1.462]./100;


%% Weights taken at every test:
mass_point = [
1 1 1 3
1 1 1 2
1 1 2 3
1 1 1 4
1 1 2 4
1 2 3 4
1 1 1 5
1 1 2 5
1 3 2 5
1 1 3 5
2 3 4 5
1 1 4 5
1 2 4 5
1 3 4 5
];

%%%%%%%%% calculations %%%%%%%%%

mass  = sum(Weights(mass_point), 2);

mass_error = zeros(1, 14);
X_error = zeros(1, 14) + 0.0000172;

for p = 1:14
    for i = 1:4
        mass_error(p) = mass_error(p) + (Weights_error(mass_point(p,i)))^2;
    end
end

mass_error = sqrt(mass_error);

gravity_force = -1.*GravityForce(mass);
gravity_force_error = -1.*(GravityForce(mass + mass_error') - GravityForce(mass)); 


%%%%%%%%% graphs %%%%%%%%%

figure 
hold on

f = fittype('a.*x + b','coefficients', {'a', 'b'});
final_fit = fit(X, gravity_force, f);
graph = plot(final_fit, 'b');
errorbar(X, gravity_force, gravity_force_error, gravity_force_error, X_error, X_error, 'color','blue','LineStyle','none', 'LineWidth', 2)

legend(graph,'K const fit', 'mg results')
grid on
box on
ylabel('MG(N)','FontSize',13)
xlabel('X(M)','FontSize',13)

hold off
