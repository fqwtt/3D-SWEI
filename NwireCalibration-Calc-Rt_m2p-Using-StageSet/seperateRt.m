function dataout=seperateRt(filein,fileout_R,fileout_t)
fidin=fopen(filein,'r');
fidout_R = fopen(fileout_R,'w');
fidout_t = fopen(fileout_t, 'w');
nline=0;
while ~feof(fidin) 
    tline=fgetl(fidin);
    nline=nline+1;
    if mod(nline, 2) == 1 
        fprintf(fidout_R,'%s\n',tline);
        dataout=tline;
    else 
        fprintf(fidout_t,'%s\n',tline);
        dataout=tline;
    end
end
fclose(fidin);
fclose(fidout_R);
fclose(fidout_t);