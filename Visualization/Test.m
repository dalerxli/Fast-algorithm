display('Solution of (red):');
tic 
SlideCut('C:\Users\nhantran.WIN\Desktop\Research\XSEDE\ScatFFT64MemOp\Test1\Results_p8000_c64000\scat9a\y.bin',20,'red.txt');
toc

display('Solution of (ie):');
tic 
SlideCut('C:\Users\nhantran.WIN\Desktop\Research\XSEDE\ScatFFT64MemOp\Test1\Results_p8000_c64000\scat9a\z.bin',20,'ie.txt');
toc

display('Solution of (ori):');
tic 
SlideCut('C:\Users\nhantran.WIN\Desktop\Research\XSEDE\ScatFFT64MemOp\Test1\Results_p8000_c64000\scat9a\FFTx.bin',20,'ori.txt');
toc
