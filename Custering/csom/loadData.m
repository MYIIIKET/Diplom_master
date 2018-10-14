function result = loadData(names)
result = struct;
for i=1:length(names)
    result.(names{i}) = load(names{i});
    result.(names{i}) = result.(names{i}).res;
end
end