function genCSV(time,cx,cy,cz,filename)

for k=1:size(cx,1)
    cx(k,:) = cx(k,:).*time.^(k-1);
    cy(k,:) = cy(k,:).*time.^(k-1);
    cz(k,:) = cz(k,:).*time.^(k-1);    
end

traj_matrix = [time; cx; cy; cz];
dlmwrite(filename,traj_matrix,'delimiter', ',' , 'precision', 10);

end