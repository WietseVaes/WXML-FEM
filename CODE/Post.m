figure(3);
subplot(131)
trisurf(elmat, x, y, u); cmocean('curl'); shading flat; title("Approximation",'Interpreter','latex','FontSize',20);
xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2)
colorbar();
subplot(132)
trisurf(elmat, x, y, usol(x,y)); cmocean('curl'); shading flat; title("Solution",'Interpreter','latex','FontSize',20);
xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2)
colorbar();
subplot(133)
trisurf(elmat, x, y, abs(usol(x,y)-u)); cmocean('curl'); shading flat; title("$|error|$",'Interpreter','latex','FontSize',20);
xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2);
colorbar();