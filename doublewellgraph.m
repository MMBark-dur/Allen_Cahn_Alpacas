syms x
F(x) = (((x^2)-1)^2)/4;
df = diff(F);
ddf = diff(diff(F));
turnpoints = solve(df == 0,x).';
minima = [];
for i = 1:length(turnpoints)
    if ddf(turnpoints(i))>0
        minima = [minima, turnpoints(i)];
    end
end
fprintf('Minima are at: u = %.1f and u = %.1f', minima(1), minima(2));
figure;
fplot(F, 'LineWidth', 3)
hold on
scatter(minima,F(minima),'or','MarkerFaceColor','r');
legend('F(u)', 'Minima' )
xlabel('u');
xlim([-2,2])
xticks(-2:1:2)
ylabel('F(u)', 'FontSize', 24);
set(gca, 'FontSize', 24)
grid on
hold off
