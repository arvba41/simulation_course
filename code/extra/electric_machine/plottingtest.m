function plottingtest(x,y,itername,xstr,ystr)

subplot(211);
plot(x,y,'LineWidth',2); grid on; hold on
legend(char(itername)); figtex(gca,1);
xlabel(char(xstr)); ylabel(char(ystr));

end
