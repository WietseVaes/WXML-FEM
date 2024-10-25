%% Compute your solution u by first building you M,P,S and b matrices and then solving the system.
[M,P,S,b] = Build(x,y,elmat,elmatbd,f,g, vx, vy);

u = (M- P + D*S)\b;