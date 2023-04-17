clear

r = linspace(0,20,1000);
d = 10;
gamma = 1;

penalty = zeros(length(r),1);
for i = 1:length(r)
    if r(i) < d
        penalty(i) = gamma*(1-(r(i)/d)^2)^3;
    else
        penalty(i) = 0;
    end
end

rho = 1;
penalty_smooth = gamma/2*(1-tanh((r-d)/rho));

penalty_smooth2 = gamma/2*(1-tanh((r-d)/1E-3));

% figure(1)
% plot(r,penalty)
% hold on
% plot(r,penalty_smooth)
% ylim([0,gamma*1.2])
% xline(d,'r--','Exclusion')
% xlabel('Distance between Spacecraft')
% ylabel('penalty')
% legend('Soft Exclusion','Absolute Exclusion')

figure(1)
plot(r,penalty_smooth,'r-')
hold on
plot(r,penalty_smooth2,'b-')
ylim([0,gamma*1.2])
xline(d,'k--','Exclusion Zone')
xlabel('Distance between Spacecraft')
ylabel('Penalty')
legend('Relaxed Exclusion','Strict Exclusion')