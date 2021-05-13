load lenses/lens.mat

conedeg=18;
cradeg=4
dAvignet = anglesVignettedLens(conedeg,cradeg,lens.exitpupil,lens.P,lens.h);
dAideal = anglesIdealLens(conedeg,cradeg);

phi=linspace(0,30,100);

for i=1:numel(phi)
   Aideal(i)=dAideal(phi(i)); 
   Avignet(i)=dAvignet(phi(i));
end
figure(5);clf;
hold on;
plot(phi,Aideal,'k','linewidth',2)
plot(phi,Avignet,'r--','linewidth',2)
