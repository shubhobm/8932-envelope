ais=importdata('aisssf.txt');
ais=ais.data;
[WX,W,f,d] = ldr(ais(:,1),ais(:,2:8),'LAD','cont','lrt','nslices',5)
[WX,W,f,d] = ldr(ais(:,1),ais(:,2:8),'LAD','cont',2,'nslices',5)


ion=importdata('ion.txt');
ion=ion.data;
[WX,W,f,d] = ldr(ion(:,1),ion(:,2:19),'LAD','disc','lrt')
[WX,W,f,d] = ldr(ion(:,1),ion(:,2:19),'LAD','disc',2)