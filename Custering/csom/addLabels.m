function dataStruct = addLabels(dataStruct, mode)
fieldNames = fieldnames(dataStruct.struct_data);
a=1;
b=0;
k=1;
label = '';
if strcmp(dataStruct.voice, 'fsp')
    label = 'F';
else
    if strcmp(dataStruct.voice, 'msp')
        label = 'M';
    end
end
for i=dataStruct.dictors
    for j=dataStruct.utterances
        b = b+length(dataStruct.struct_data.(fieldNames{k}));
        k=k+1;
    end
    d = a:b;
    if strcmp(mode, 'dictor')
        index = num2str(i);
    else
        if strcmp(mode, 'data')
            index = '';
        end
    end
    dataStruct.som_data = som_label(dataStruct.som_data,...
        'add', d, strcat(label,index));
    a = b+1;
end
end