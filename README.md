Řešení LODR1 pomocí Eulerovy metody


% Parametry
t0 = 0;        % Počáteční čas
x0 = 1;        % Počáteční podmínka
h = 1;       % Krok
t_end = 10;     % Konec intervalu
t = t0:h:t_end; % Vektor časových okamžiků

odeFcn = @(t,x) -0.5*x;

Analytické řešení
x_analytic = exp(-0.5 * t);

Eulerova metoda
x_euler = zeros(size(t));
x_euler(1) = x0;
for i = 1:length(t)-1
    a_euler = odeFcn(t, x_euler(i));%-0.5 * x_euler(i)
    x_euler(i+1) = x_euler(i) + h * a_euler;
end

Lichoběžníková metoda
x_lichobeznikova = zeros(size(t));
x_lichobeznikova(1) = x0;
for i = 1:length(t)-1
    f1 = odeFcn(t, x_lichobeznikova(i));
    f2 = -0.5 * (x_lichobeznikova(i) + h * f1);
    x_lichobeznikova(i+1) = x_lichobeznikova(i) + h * (f1 + f2) / 2;
end

Runge-Kutta 4. řádu (ode4)
x_tk = zeros(size(t));
x_tk(1) = x0;
for i = 1:length(t)-1
    a1 = -0.5 * x_tk(i);
    a2 = -0.5 * (x_tk(i) + h * a1 / 2);
    a3 = -0.5 * (x_tk(i) + h * a2 / 2);
    a4 = -0.5 * (x_tk(i) + h * a3);
    x_tk(i+1) = x_tk(i) + (h / 6) * (a1 + 2*a2 + 2*a3 + a4);
end

Výpočet MSE
mse_euler = mean((x_euler - x_analytic).^2);
mse_trap = mean((x_lichobeznikova - x_analytic).^2);
mse_rk = mean((x_tk - x_analytic).^2);

Vykreslení výsledků
figure;
hold on;
plot(t, x_analytic, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytické řešení', 'Color','yellow');
plot(t, x_euler, 'r--o', 'DisplayName', 'Eulerova metoda');
plot(t, x_lichobeznikova, 'g-.s', 'DisplayName', 'Lichoběžníková metoda');
plot(t, x_tk, 'b-x', 'DisplayName', 'Runge-Kutta metoda');
xlabel('Čas (t)');
ylabel('x(t)');
title('Řešení diferenciální rovnice pomocí numerických metod');
legend;
grid on;
hold off;

% Zobrazení MSE
fprintf('MSE Eulerova metoda: %.6f\n', mse_euler);
fprintf('MSE Lichoběžníková metoda: %.6f\n', mse_trap);
fprintf('MSE Runge-Kutta metoda: %.6f\n', mse_rk);
