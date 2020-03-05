AreaLength=1000;
SensorNum=400;
SensorX=rand(1,40)*1000;
SensorY=rand(1,40)*1000;
figure(1);
scatter(SensorX,SensorY,8);
for i=1:40
    Str = int2str(i);
    text(SensorX(i)+10,SensorY(i)+10,Str,'FontName','Times New Roman','FontSize',12);
end