from matplotlib import pyplot as plt
import re


if __name__ == '__main__':
    print("\n ---- Plotting forces ----\n")

    # File name
    filename='postProcessing/forces1/0/forces.dat'


    # Set empty lists
    t = []
    fpx = []; fpy = []; fpz = []
    fvx = []; fvy = []; fvz = []

    # Read file line by line and populate arrays
    pipefile=open(filename,'r')
    lines = pipefile.readlines()
    for line in lines:
        line=re.sub('\(','',line)
        line=re.sub('\)','',line)
        ls=line.split()
        if ls[0] != '#':
                t.append(float(ls[0]))
                fpx.append(float(ls[1]))
                fpy.append(float(ls[2]))
                fpz.append(float(ls[3]))
                fvx.append(float(ls[4]))
                fvy.append(float(ls[5]))
                fvz.append(float(ls[6]))


    # Clean up data files
    # Create dictionary of values
    # And add in latest data
    tdict={}
    fx={}
    fy={}
    fz={}
    for i in range(len(t)):
        tdict[t[i]]=t[i]
        fx[t[i]] = fpx[i] + fvx[i]
        fy[t[i]] = fpy[i] + fvy[i]
        fz[t[i]] = fpz[i] + fvz[i]


    time=list(fx.keys())
    forceX=list(fx.values())
    forceY = list(fy.values())


    plt.xlabel('Time (s)')
    plt.ylabel('Force, X (N)')
    plt.grid(True)
    plt.plot(time,forceX,'-')
    plt.savefig('forceX.pdf',bbox_inches='tight')
    plt.close()

    plt.xlabel('Time (s)')
    plt.ylabel('Force, Y (N)')
    plt.grid(True)
    plt.plot(time,forceY,'-')
    plt.savefig('forceY.pdf',bbox_inches='tight')
    plt.close()


