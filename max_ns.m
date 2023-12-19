function max_nodal_stress = max_ns(Stress,column)

[rows,columns] = size(Stress);

max_value = 0;
row_value = 1;
for i = 1:rows
      if abs(Stress(i,column)) > abs(max_value)
          max_value = Stress(i,column);
          row_value = i;
      end
end

max_nodal_stress = [row_value,max_value];


end