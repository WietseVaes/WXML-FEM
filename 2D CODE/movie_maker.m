function [] = movie_maker(u,x,y, elmat,t, sp, Gif_name,s)

maxvals = max(max(u))
minvals = min(min(u))

% Set figure size to fixed dimensions (1120x840 pixels)
f = figure(365); hold on
f.Position = [100, 100, 1120, 840];
f.Color = [1 1 1]; 
for i1 = sp:length(t)
    uplot = u(:,i1);

    trisurf(elmat,x,y,uplot); shading interp;cmocean('topo'); view(2);
    title(['Solution at time $t =~$',num2str(t(i1),'%.3f')],'Interpreter','latex','FontSize',20)
    xlabel('$x$','Interpreter','latex','FontSize',20); ylabel('$y$','Interpreter','latex','FontSize',20)
    zlim([minvals,maxvals]); 
    colorbar();
    caxis([minvals,maxvals]); % Makes the colors not jump
    symabsval = max(abs(minvals),abs(maxvals));
    zlim([-symabsval,symabsval]); 
    colorbar();
    caxis([-symabsval,symabsval]); % Makes the colors not jump
    axis square 
    grid off
    axis off;
    hold on; 
    plot(s.Z,'linewidth',2,'k')
    Image = getframe(gcf);
    
    axis square;

    outputDir = 'Results';

    outputFile = fullfile(outputDir, ['im_', num2str(i1 - sp), '.png']);
    imwrite(Image.cdata, outputFile, 'png');
    clf;
end
Gif_writing(Gif_name);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%subfunctions%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GIF maker

function [] = Gif_writing(gifFile)
outputDir = 'Results';
outgifFile = fullfile(outputDir,gifFile);
i = 1;
while isfile(fullfile(outputDir, ['im_', num2str(i - 1), '.png']))
    images{i} = imread(fullfile(outputDir, ['im_', num2str(i - 1), '.png']));
    delete(fullfile(outputDir, ['im_', num2str(i - 1), '.png']))
    i = i + 1;
end

n = i-1;
writerObj = VideoWriter(outgifFile,'MPEG-4');
writerObj.FrameRate = 20;
secsPerImage = ones(n) ;
open(writerObj);

for u=1:n
    % convert the image to a frame
    frame = im2frame(images{u});
    for v=1:secsPerImage(u)
        writeVideo(writerObj, frame);
    end
end

close(writerObj);
end