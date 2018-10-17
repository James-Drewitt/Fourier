###########################################
#  ___             _               _   _  #
# | __|__ _  _ _ _(_)___ _ _  __ _/ | / | #
# | _/ _ \ || | '_| / -_) '_| \ V / |_| | #
# |_|\___/\_,_|_| |_\___|_|    \_/|_(_)_| #  
#                                         #
# Fourier v1.1                            #
# 29/04/2013                              #
# James Drewitt (james.drewitt@gmail.com) #
###########################################
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import numpy as np
from Tkinter import *
#################
###Subroutines###
#################
def FFT(SQ, rho, M, maxq):
    print '*** Fast Fourier Transform ***'
    nrows=SQ.shape[0]
    ntrans=(2**M)-1
    nhalf=(ntrans+1)/2
    qstep=SQ[1,0]-SQ[0,0]
    qmax=nhalf*qstep
    nmax=int(np.math.ceil(maxq/qstep))
    rstep=np.math.pi/qmax
    F=np.zeros((ntrans+1,1))
    k=1
    ncut=min(nmax,nrows)
    while k<ncut:
        F[k]=SQ[k,1]*k*qstep
        F[ntrans-k+1]=-F[k]
        k=k+1
    FT=np.zeros((len(F),1))
    FT=np.imag(np.fft.fft(F, axis=0)/ntrans)
    r=np.zeros((nhalf,))
    G=np.zeros((nhalf,2))
    for m in range(nhalf):
        r[m]=m*rstep
        G[m,0]=r[m]
        G[m,1]=-FT[m]*(qmax/(2*np.math.pi*np.math.pi*float(rho)*G[m,0]))
    return G
def BFFT(Gr, rho, M):
    print '*** Back FFT ***'
    nrows=Gr.shape[0]
    ntrans=(2**M)-1
    nhalf=(ntrans+1)/2
    rstep=Gr[1,0]-Gr[0,0]
    rmax=Gr[nrows-1,0]+rstep
    qstep=np.math.pi/rmax
    F=np.zeros((ntrans+1,1))
    k=1
    while k<nrows:
        F[k]=Gr[k,1]*k*rstep
        F[ntrans-k+1]=-F[k]
        k=k+1
    FT=np.imag(np.fft.fft(F, axis=0))/ntrans
    q=np.zeros((nhalf,1))
    BT=np.zeros((nhalf,2))
    for m in range(1,nhalf):
        q[m]=(m)*qstep
        BT[m,0]=q[m]
        BT[m,1]=-FT[m]*4*np.math.pi*float(rho)*rmax/q[m]
    return BT
####################
# INPUT PARAMETERS #
####################
root=Tk()
root.title("Fourier v1.1")
root.minsize(250,100);root.resizable(0,0)
def letsgo():
    global name,path,rho
    name=inp1.get();path=inp2.get();rho=inp3.get()
    root.quit();root.destroy()
l1=Label(root,text="File name: ").grid(row=0,column=1, sticky=W)
l2=Label(root,text="File directory: ").grid(row=1,column=1,sticky=W)
l3=Label(root,text="Number density: ").grid(row=2,column=1,sticky=W)
inp1=Entry(root);inp2=Entry(root);inp3=Entry(root);inp4=Entry(root)
inp1.insert(0,"CA_2000_SQ.dat")
inp2.insert(0,"O://00_Python//FT//")
inp3.insert(0,"0.074")
inp1.grid(row=0,column=2,sticky=W),inp2.grid(row=1,column=2,sticky=W)
inp3.grid(row=2,column=2,sticky=W)
w=Button(root, text="FFT", width=10,command=letsgo,foreground="blue")
w.grid(row=4,column=1,columnspan=2)
root.mainloop()
#
filepath=path+name
data=np.loadtxt(filepath)
print data
###########################
# PLOT DATA AND TRANSFORM #
###########################
def _quit():
    root.quit();root.destroy()
def foo():
    m1.destroy();m2.destroy()
    bar()
def bar():
    global m1, m2
    m1=Frame(midframe);m1.pack(side=LEFT)
    m2=Frame(midframe);m2.pack(side=LEFT)
    lowq=int(e0.get());rc1=float(e2.get())
    rc2=float(e3.get());rho=float(e1.get())
    cospar=int(e4.get());SRdispl=float(e7.get())
    sc=float(esc.get())
    qstep=0.05
    print data
    SQ2=data
    print SQ2
    Qmax=SQ2[SQ2.shape[0]-1,0];print "Qmax=";print Qmax
    SQ2[0:lowq,1]=SQ2[lowq,1]
    SQ2[:,1]=SQ2[:,1]*sc
    G=FFT(SQ2, rho, 12, Qmax)
    G1=np.zeros((rc1,2))
    G1[:rc1,0]=G[:rc1,0];G1[:rc1,1]=G[:rc1,1];
    G[0:rc1,1]=-1
    BT1=BFFT(G, rho, 12)
    BT1[:,1]=BT1[:,1]+SRdispl
    BT1=np.delete(BT1,np.s_[SQ2.shape[0]-1:],axis=0)
    OPBT1=path+name+"-sq.dat"
    np.savetxt(OPBT1,BT1)
    NOCOS=np.zeros((BT1.shape[0],2))
    NOCOS[:,0]=BT1[:,0]
    NOCOS[:,1]=BT1[:,1]
    if cospar>0:
        cosmax=BT1.shape[0];print cosmax
        cosmin=cosmax-cospar
        for K in range(cosmin,cosmax):
            J=K-cosmin
            BT1[K,1]=BT1[K,1]*0.5*(1.0+np.math.cos(J*np.math.pi/(cospar-1)))
    G2=FFT(BT1, rho, 12, Qmax)
    OPGR2=path+name+"-Gr.dat"
    np.savetxt(OPGR2,G2)
    G3=np.zeros((rc2,2))
    G3[:rc2,0]=G2[:rc2,0];G3[:rc2,1]=G2[:rc2,1]
    G2[0:rc2,1]=-1
    BT2=BFFT(G2,rho,12)
    BT2=np.delete(BT2,np.s_[SQ2.shape[0]-1:],axis=0)
    OPBT2=path+name+"-BT.dat"
    np.savetxt(OPBT2,BT2)
    f=plt.figure(dpi=100)
    plt.plot(SQ2[:,0],SQ2[:,1],'c--',NOCOS[:,0],NOCOS[:,1],'k-',BT2[:,0],BT2[:,1],'r-')
    plt.xlim( 0, Qmax )
    global canvas
    canvas=FigureCanvasTkAgg(f, m1)
    toolbar=NavigationToolbar2TkAgg(canvas,m1) 
    toolbar.update()
    canvas.get_tk_widget()
    canvas._tkcanvas.pack(side=LEFT)
    canvas._tkcanvas.configure(width=500,height=500)
    f2=plt.figure(dpi=100)
    global canvas2
    canvas2=FigureCanvasTkAgg(f2, m2)
    toolbar2=NavigationToolbar2TkAgg(canvas2,m2) 
    toolbar2.update()
    canvas2.get_tk_widget()
    canvas2._tkcanvas.pack(side=LEFT)
    canvas2._tkcanvas.configure(width=500,height=500)
    plt.axhline(y=-1,color='b')
    plt.plot(G2[:,0],G2[:,1],'k-')
    plt.autoscale(True,'y',False)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,y1-np.abs(y1)*0.2,y2+np.abs(y2)*0.05))
    plt.autoscale(False)
    plt.plot(G1[:,0],G1[:,1],'c--',G3[:,0],G3[:,1],'r--')
    plt.xlim(0, 8)
    canvas.draw();canvas2.draw();
    SR=0
    print len(BT1)
    for k in range(len(BT1)):
        SR=SR+(BT1[k,1]*BT1[k,0]*BT1[k,0]*(BT1[2,0]-BT1[1,0]))  
    e5.delete(0,END);e5.insert(0,SR)
    SRtheo=-2*np.math.pi*float(rho)
    e6.delete(0,END);e6.insert(0,SRtheo)
root=Tk();root.resizable(0,0)
frame=Frame(root);frame.pack(side=TOP)
midframe=Frame(root);midframe.pack(side=TOP)
botframe=Frame(root);botframe.pack(side=BOTTOM)
l0=Label(frame,text="low-Q cut");l0.pack(side=LEFT)
e0=Entry(frame);e0.insert(0,"3");e0.pack(side=LEFT)
l1=Label(frame,text="Number density: ");l1.pack(side=LEFT)
e1=Entry(frame);e1.insert(0,rho);e1.pack(side=LEFT)
B1=Button(frame,text="FFT",command=foo, foreground="blue");B1.pack(side=LEFT)
B2=Button(frame,text="Quit",command=_quit, foreground="red");B2.pack(side=RIGHT)
#
lsc=Label(botframe,text="Scale factor: ");lsc.pack(side=LEFT)
esc=Entry(botframe,width=5);esc.insert(0,"1");esc.pack(side=LEFT)
#
l2=Label(botframe,text="Slope correction: ");l2.pack(side=LEFT)
e2=Entry(botframe,width=5);e2.insert(0,"12");e2.pack(side=LEFT)
l3=Label(botframe,text="r cut:");l3.pack(side=LEFT)
e3=Entry(botframe,width=5);e3.insert(0,"22");e3.pack(side=LEFT)
l4=Label(botframe,text="Cosine smooth factor:");l4.pack(side=LEFT)
e4=Entry(botframe,width=5);e4.insert(0,"0");e4.pack(side=LEFT)
l5=Label(botframe,text="Sum rule:");l5.pack(side=LEFT)
e5=Entry(botframe,width=10);e5.insert(0,"0");e5.pack(side=LEFT)
l6=Label(botframe,text="-2*pi*n0:");l6.pack(side=LEFT)
e6=Entry(botframe,width=10);e6.insert(0,"0");e6.pack(side=LEFT)
l7=Label(botframe,text="displacement:");l7.pack(side=LEFT)
e7=Entry(botframe,width=10);e7.insert(0,"0");e7.pack(side=LEFT)
bar()
###############
root.mainloop()
###############
plt.close('all')
