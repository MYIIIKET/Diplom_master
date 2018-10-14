function dataStruct = loadDataStruct(dictors, utterances,...
    variance, data, noise, voice, noiseValue, mode)
if strcmp(noise, '_orig')
    noiseValue = '';
end
names = getFileNames(dictors,...
    utterances, voice, variance, data, noise, noiseValue);
struct_data = loadData(names);

flat_data = [];
for i=1:length(names)
    flat_data = [flat_data, struct_data.(names{i})];
end

som_data = som_data_struct(flat_data');
som_data = som_normalize(som_data, 'var');

dataStruct = struct;
dataStruct.dictors = dictors;
dataStruct.utterances = utterances;
dataStruct.variance = variance;
dataStruct.data = data;
dataStruct.voice = voice;
dataStruct.noise = noise;
dataStruct.noiseValue = noiseValue;
dataStruct.names = names;
dataStruct.struct_data = struct_data;
dataStruct.flat_data = flat_data;
dataStruct.som_data = som_data;
dataStruct = addLabels(dataStruct, mode);
end