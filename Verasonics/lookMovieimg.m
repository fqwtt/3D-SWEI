% obj = VideoReader('G:\WT\ºÙ«–≤®\0615\◊Û°™”“Àƒ¥ŒºÙ«–≤®.mp4');% ‰»Î ”∆µŒª÷√
% frame_num=obj.NumFrames;
% data=zeros(371-35+1,653-96+1,frame_num);
% for i=1:frame_num
%     frame=read(obj,i); %#ok<*VIDREAD>
%     frame=frame(35:371,96:653);
%     data(:,:,i)=frame;
% end
movie=dis;
[~,~,frames]=size(movie);
figure
for i=1:1:frames
    
%     colormap('gray')
    
    imagesc((movie(:,:,i)));
    drawnow
    title('Displacement');
    xlabel('lateral position(wavelength)');ylabel('axial position(wavelength)');
    set(gca,'XTick',linspace(1,200,11),'XTickLabel',round(linspace(-50,50,11)),...
            'YTick',linspace(1,400,7),'YTickLabel',round(linspace(0,100,7)))
    

end
% for i=100:150
%     plot(-squeeze(dis(200,i,1:80)));
% end
% cmap = parula(256);
% sliceViewer(dis,'Colormap',cmap,'DisplayRange',[-10,10],'ScaleFactors',[3,2,1])