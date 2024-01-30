clc, clear, close all;

n = [1;2;3;4;5];
A = [42; 162; 642; 2562; 10242];
Harmonic_Eqn_1 = [0.8452;0.2297;0.0586;0.0147;0.0037];
Harmonic_Eqn_2 = [2.2936;0.6766;0.1761;0.0445;0.0112];

err_l2 = table(n, Harmonic_Eqn_1, Harmonic_Eqn_2)
set(gca,'linewidth',3,'fontsize',25,'TickLabelInterpreter','latex')

subplot(2,1,1);
loglog(A, Harmonic_Eqn_1, '-*')
hold on
p = polyfit(log(A),log(Harmonic_Eqn_1),1);
z = polyval(p,log(A));
hold on
loglog(A,exp(z), "r--")
hold on
xlabel('Precision of Triangular Mesh','interpreter','latex','fontsize',30)
ylabel('Error','interpreter','latex','fontsize',30)
title('Loglog plot for error (Harmonic Equation 1)','interpreter','latex','fontsize',35)
legend('Loglog plot','Best fit curve','interpreter','latex','fontsize',35) 

subplot(2,1,2);
loglog(A, Harmonic_Eqn_2, '-*')
hold on
p2 = polyfit(log(A),log(Harmonic_Eqn_2),1);
z2 = polyval(p2,log(A));
hold on
loglog(A,exp(z2), "r--")
xlabel('Precision of Triangular Mesh','interpreter','latex','fontsize',30)
ylabel('Error','interpreter','latex','fontsize',30)
title('Loglog plot for error (Harmonic Equation 2)','interpreter','latex','fontsize',35)
legend('Loglog plot','Best fit curve','interpreter','latex','fontsize',35)