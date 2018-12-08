


clear;

for i=1:100
   [INL(i,:), DNL(i,:)] = getDNLINL();

end

sigma_DNL = std(DNL);
sigma_INL = std(INL);

sigma_3_DNL = 3.*max(sigma_DNL);
sigma_3_INL = 3.*max(sigma_INL);


figure(1)
subplot(2,1,1);
plot(sigma_DNL,'DisplayName','DNL','LineWidth',2);
ylabel('\sigma_{DNL} [LSB]','FontSize',12,'FontWeight','bold');
xlabel('Output Code','FontSize',12,'FontWeight','bold');

grid on
legend('show');
xlim([0,2^11]);

subplot(2,1,2);
plot(sigma_INL,'DisplayName','INL','LineWidth',2);
ylabel('\sigma_{INL} [LSB]','FontSize',12,'FontWeight','bold');
xlabel('Output Code','FontSize',12,'FontWeight','bold');

grid on
legend('show');
xlim([0,2^11]);
