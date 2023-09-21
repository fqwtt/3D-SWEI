function dis = lowpass2(dis,thread)
    [depth,width,~]=size(dis);
    mul=zeros(depth,width);
    for i = 1 : depth
        for j = 1 : width
            if sqrt((i-(1+depth)/2)^2+(j-(1+width)/2)^2) < thread
                mul(i,j)=1;
            end
        end
    end
    spectrum_dis = fftshift(fft2(dis));
    spectrum_dis = spectrum_dis.* mul;
    dis=real(ifft2(ifftshift(spectrum_dis)));
end

