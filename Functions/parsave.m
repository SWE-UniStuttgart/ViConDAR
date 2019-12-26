% Cheating parfor to convince it to save...

function parsave(tupe,path,variable)
windfield=variable;
save(tupe,path,'windfield')
end