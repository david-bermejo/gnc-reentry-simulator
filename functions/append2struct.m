function res = append2struct(structObj, fields, values, dim)
    res = structObj;
    for i=1:dim
        res.(fields{i}) = values{i};
    end
end