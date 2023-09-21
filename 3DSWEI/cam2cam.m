cs2c1=zeros(3,3,8);
cs2c1(:,:,1)=eye(3);
cs2c1(:,:,2)=cs2c1(:,:,1)*(Data.w2c(:,:,1,1)*Data.Rw2Lw/Data.w2c(:,:,2,1));
cs2c1(:,:,3)=cs2c1(:,:,2)*(Data.w2c(:,:,1,2)/Data.Rw2Lw/Data.w2c(:,:,2,2));
cs2c1(:,:,4)=cs2c1(:,:,3)*(Data.w2c(:,:,1,3)/Data.Rw2Lw/Data.w2c(:,:,2,3));
cs2c1(:,:,5)=cs2c1(:,:,4)*(Data.w2c(:,:,1,4)*Data.Rw2Lw/Data.w2c(:,:,2,4));
cs2c1(:,:,6)=cs2c1(:,:,5)*(Data.w2c(:,:,1,5)/Data.Rw2Lw/Data.w2c(:,:,2,5));
cs2c1(:,:,7)=cs2c1(:,:,6)*(Data.w2c(:,:,1,6)/Data.Rw2Lw/Data.w2c(:,:,2,6));
cs2c1(:,:,8)=cs2c1(:,:,7)*(Data.w2c(:,:,1,7)*Data.Rw2Lw/Data.w2c(:,:,2,7));

fid=fopen('cs2c1.txt','w');
for k=1:8
    for i=1:3
        for j=1:3
            fprintf(fid,'%.6g ',cs2c1(i,j,k));
        end
    end
    fprintf(fid,'\n');
end

