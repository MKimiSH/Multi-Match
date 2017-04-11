function ShowInstance(I1,I2,prefixName,templateMask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('templateMask','var'))
    templateMask = ones(size(I1,1), size(I1,2));
end

I1(templateMask==0) = 0;

[h1,w1,d] = size(I1);
[h2,w2,d] = size(I2);

r1x = 0.5*(w1-1);
r1y = 0.5*(h1-1);
r2x = 0.5*(w2-1);
r2y = 0.5*(h2-1);

[xs,ys] = meshgrid(1:w1,1:h1);

% take the complete set of pixels (not a sample)
xs = xs(:)'; 
ys = ys(:)';





    
% a mask on ones
%I3 = I1;
%I3(:) = 1;
I3 = templateMask;
% transformed mask



%% summarizing figure (image)
tic
fullscreen = get(0,'ScreenSize');
figure()
set(gcf,'Position',[0.15*fullscreen(3) 0.1*fullscreen(4) 0.8*fullscreen(3) 0.75*fullscreen(4)]);
% GL.hh = hh;
set(gcf,'color','w');
set(gcf,'name',[prefixName, ': FM result']);

%I1
subplot 221; hold off;
% I1padded = padarray(I1,[max(0,h2-h1),max(0,floor(0.5*(w2-w1)))],1,'post'); 
% I1padded = padarray(I1padded,[0,max(0,ceil(0.5*(w2-w1)))],1,'pre'); 
white_I1 = I1;
white_I1(~templateMask) = 1;
I1padded_the_size_of_I2 = padarray(white_I1,[round((h2-h1)/2),round((w2-w1)/2)],1);
imshow(I1padded_the_size_of_I2); hold on;
title(['TemplateSize: ' num2str(h1) 'x' num2str(w1) ...
        '.   TV: ' num2str(GetTotalVariation(I1),'%.3f') '.  STD: ' num2str(std(I1(:)),'%.3f')]);

% I2
subplot 222; hold off;
imshow(I2); hold on;
return

% secondSTR = sprintf('\nwarped template (magenta)');
% if (~isempty(optA))
%     plot([cornerOPTxs cornerOPTxs(1)],[cornerOPTys cornerOPTys(1)], '*-g');
%     secondSTR = sprintf('%s, ground truth (green)',secondSTR);
% end
% plot([cornerAxs cornerAxs(1)],[cornerAys cornerAys(1)], '*-m');
% % title(['TargetSize: ' num2str(h2) 'x' num2str(w2) '.     outside: ' num2str((numPoints-length(insideInds))/numPoints,'%.3f') ' L1-dist: ' num2str(fullError,'%.4f') ' overlapErr: ' num2str(overlapError,'%.4f')]);
% title(['TargetSize: ' num2str(h2) 'x' num2str(w2) '.   L1-dist: ' num2str(fullError,'%.3f') ...
%     '   overlapErr: ' num2str(overlapError,'%.3f') secondSTR]);
% 
return
% 
% % query result
% subplot 223; hold off;
% maskedtranstarget = ~transI3mask.*ones(size(transI3mask)) + transI3mask.*transtarget;
% imshow(maskedtranstarget(yLimits(1):yLimits(end),xLimits(1):xLimits(end),:)); hold on;
% title('the warped template');
% 
% % ground truth
% subplot 224; hold off;
% transI3mask = double(transI3mask);
% target =  ~transI3mask.*ones(size(transI3mask)) + transI3mask.*I2;
% 
% imshow(target(yLimits(1):yLimits(end),xLimits(1):xLimits(end))); hold on;
% title('background target');
% 
% 
% return
% 
% 
% 
% 
% 
% 
% 
