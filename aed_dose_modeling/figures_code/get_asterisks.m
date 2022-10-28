function ast = get_asterisks(p,adjust)

if p < 0.001/adjust
    ast = '***';
elseif p < 0.01/adjust
    ast = '**';
elseif p < 0.05/adjust
    ast = '*';
else
    ast = 'ns';
end
    

end