function m = RingDataset(  ri,ro,c,samples_num)
% Создаем набор данных
m = rand(2,samples_num)*2*ro-ro;
c = repmat(c,1,samples_num);
ee = m+c;
een = sqrt(ee(1,:).^2+ee(2,:).^2);
ind = find( (een<ro) & (een>ri) );
m = m( :, ind );

while(size(m,2)<samples_num)
   mt = rand(2,1)*2*ro-ro; 
   ee = mt+c(:,1);
   een = sqrt(ee(1,:).^2+ee(2,:).^2);
   if ((een<ro) & (een>ri) )
       m = [m mt];
   end;
end;

end