
% figure(3);
% subplot(131)
% trisurf(elmat, x, y, u); cmocean('curl'); shading flat; title("Approximation",'Interpreter','latex','FontSize',20);
% xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2)
% colorbar();
% subplot(132)
% trisurf(elmat, x, y, usol(x,y)); cmocean('curl'); shading flat; title("Solution",'Interpreter','latex','FontSize',20);
% xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2)
% colorbar();
% subplot(133)
% trisurf(elmat, x, y, abs(usol(x,y)-u)); cmocean('curl'); shading flat; title("$|error|$",'Interpreter','latex','FontSize',20);
% xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2);
% colorbar();


Gif_name = 'Approx3';
sp = 1 ;

movie_maker(u,x,y, elmat, t, sp, Gif_name)

% Gif_name = 'Solution';
% movie_maker(usol,x,y, elmat, t, sp, Gif_name)
% 
% Gif_name = 'error';
% movie_maker(abs(u-usol),x,y, elmat, t, sp, Gif_name)

% figure(3);
% subplot(131)
% trisurf(elmat, x, y, u(:,2)); cmocean('curl'); shading flat; title("Approximation",'Interpreter','latex','FontSize',20);
% xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2)
% colorbar();
% subplot(132)
% trisurf(elmat, x, y, usol(:,2)); cmocean('curl'); shading flat; title("Solution",'Interpreter','latex','FontSize',20);
% xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2)
% colorbar();
% subplot(133)
% trisurf(elmat, x, y, abs(usol(:,2)-u(:,2))); cmocean('curl'); shading flat; title("$|error|$",'Interpreter','latex','FontSize',20);
% xlabel("$x$",'Interpreter','latex','FontSize',20); ylabel("$y$",'Interpreter','latex','FontSize',20); view(2);
% colorbar();
