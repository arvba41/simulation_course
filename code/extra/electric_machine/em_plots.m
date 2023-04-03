function em_plots(t,y,Te,wout,vs)

subplot(231);
plot(t,y(1,:)); grid on; ylabel('$i_{d}(t)$ [A]'); xlabel('$t$');
ax(1) = figtex(gca); hold on; ylim([0 6000])
subplot(232); 
plot(t,y(2,:)); grid on; ylabel('$i_{q}(t)$ [A]'); xlabel('$t$');
ax(2) = figtex(gca); hold on; ylim([0 15000])

subplot(233); 
plot(t,abs(y(1,:) + 1j*y(2,:))); grid on; ylabel('$\bar i_{s}(t)$ [A]'); 
xlabel('$t$'); ax(3) = figtex(gca); hold on; ylim([0 15000])

subplot(234); 
plot(t,wout); grid on; ylabel('$\omega_{r(e)}(t)$ [rpm]');
xlabel('$t$'); ax(4) = figtex(gca); hold on; ylim([0 1500])

subplot(235); 
plot(t,Te); grid on; ylabel('$T_e(t)$ [Nm]'); xlabel('$t$'); 
ax(5) = figtex(gca); hold on; ylim([0 15000])

subplot(236); 
plot(t,abs(vs)); grid on; ylabel('$v_{dq}(t)$'); xlabel('$t$');
ax(6) = figtex(gca); hold on; 

linkaxes(ax,'x')