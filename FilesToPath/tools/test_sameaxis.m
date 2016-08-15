  subplot(3,1,1); plot(randn(100,1),randn(100,1),'x'); ylabel('QF')
  subplot(3,1,2); plot([-1 0 .5 1],[0 100 100 10],'x-'); ylabel('HT');
  subplot(3,1,3); plot(randn(100,1),randn(100,1)*33,'x'); ylabel('DV');
  samexaxis('abc','xmt','on','ytac','join','yld',1)