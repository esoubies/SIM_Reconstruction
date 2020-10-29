function setConfigFilePSF(Nx,Ny,Nz,res_xy,res_z,Na,lamb,configFile)

% Read txt into cell A
fid = fopen(configFile);
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

% Change cell A
A{3} = sprintf('Lambda=%f',lamb);
A{5} = sprintf('NA=%f',Na);
A{6} = sprintf('NX=%f',Nx);
A{7} = sprintf('NY=%f',Ny);
A{8} = sprintf('NZ=%f',Nz);
A{9} = sprintf('ResAxial=%f',res_z);
A{10} = sprintf('ResLateral=%f',res_xy);

% Write cell A into txt
fid = fopen(configFile, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
end
