import numpy as np
import PIL
import matplotlib.pyplot as plt
import math 
from mpl_toolkits import mplot3d
import matplotlib.animation as animation
fig = plt.figure(figsize = (20,20))
ax = plt.axes(projection='3d')
ax.grid()
cos = np.cos
sin = np.sin
arcsin = np.arcsin
arccos = np.arccos
f = 0
e = 0.08874
i = 7.14043 #in degrees
a = 2.36179
w = 151.19853
om = 103.85136 
pi = math.pi
ims = []

def calcAngle(sinx, cosx):
    if sinx > 0 and  cosx > 0:
        return arcsin(sinx)
    elif sinx > 0 and cosx < 0:
        return arccos(cosx)
    elif sinx < 0 and cosx < 0:
        return pi + arcsin(-1 * sinx)
    else: return arcsin(sinx)


def prepareOrbit(f, e, i, a, w, om, col):
    i = math.radians(i)
    w = math.radians(w)
    om = math.radians(om)
    p = a* (1 - e**2)
    tan = np.tan
    sini = sin(i)
    cosi = cos(i)
    inc = 2 * math.pi / 10000
    x, y, z = [list(range(10000)) for i in range(3)]
    
    for i in range(10000):
        r = p/(1 + e *cos(f))
        sinb = sini * math.sin(w + f)
        tanlom = cosi * tan(w + f)
        lom = np.arctan(tanlom)
        if (sin(lom)* sin(w+f)) <  0:
            lom += pi
        l = lom + om
        b = np.arcsin(sinb)
    
        x[i] = r * cos(b)*sin(l)
        y[i] = r * cos(b)*cos(l)
        z[i] = r * sin(b)
        f += inc
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    ax.plot(x, y, z, color=col)
    
def calcE(M, e):
    E = M
    sinE = sin(E)
    M = M%(2 * pi)
    #print(M)
    for i in range(100):
        E = M + e * sin(E) 
 #    print(f'E:{E} ')
    return E

def updateRF(a, e, t, tau, Per):
    Per = Per * 365.25
    M = 2 * pi/Per * (t - tau)
   # print(f'M: {M}, Per: {Per}, t:{t}')
    E_eSinE = M
    E = calcE(M, e)
    rnow = a*(1 - e *cos(E))
    cosf = (cos(E) - e)/(1 - e * cos(E))
    sinf = math.sqrt(1 - e**2) * sin(E)/(1 - e * cos(E))
    f = calcAngle(sinf, cosf) 
    return (rnow, f)


def precalcPos(e, i, a, w, om, dt, tau):
    i = math.radians(i)
    w = math.radians(w)
    om = math.radians(om)
    p = a* (1 - e**2)
    tan = np.tan
    Per = math.sqrt(a**3)
    sini = sin(i)
    cosi = cos(i)
    t = 0
    r, f = updateRF(a, e, t, tau, Per) 
    x, y, z = [list(range(10000)) for i in range(3)]
    for i in range(10000):
        r, f = updateRF(a,e, t, tau, Per)
        sinb = sini * math.sin(w + f)
        tanlom = cosi * tan(w + f)
        lom = np.arctan(tanlom)
        if (sin(lom)* sin(w+f)) <  0:
            lom += pi
        l = lom + om
        b = np.arcsin(sinb)
    
        x[i] = r * cos(b)*sin(l)
        y[i] = r * cos(b)*cos(l)
        z[i] = r * sin(b)
        t += dt
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    return (x, y, z)
    
#bm = PIL.Image.open('backgr.jpeg')
#bm = np.array(bm.resize([int(d/5) for d in bm.size]))/256.
#lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
#lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 
#x = np.outer(np.cos(lons), np.cos(lats)).T
#y = np.outer(np.sin(lons), np.cos(lats)).T
#z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
#ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors = bm)
#x, y, z = 10*(np.random.rand(3,100)-0.5)
#ax.scatter(x, y, z, s=0.1, c='w')
#fig.set_facecolor('black')
#ax.set_facecolor('black')
#def init(orbs):
#    for i in range(len(orbs)): 

#lines = [plt.plot([], [], [])[0] for i in range(len(params))]
redDot, = plt.plot([0], [0], [0], 'ro')
orDot, = plt.plot([0], [0], [0], 'go')
blDot, = plt.plot([0], [0], [0], 'bo')
bkDot, = plt.plot([0], [0], [0], 'ko')
#print(redDot)
def animate(i):
    ax.set_title(f'{i*10} days')
    orDot.set_data_3d(x1[i], y1[i], z1[i])
    redDot.set_data_3d(x3[i], y3[i], z3[i])
    blDot.set_data_3d(x2[i], y2[i], z2[i])
    bkDot.set_data_3d(x4[i], y4[i], z4[i])
    return blDot, orDot, redDot, bkDot


def animate2(i):
    for line in lines:
        line.set_data_3d(*coords, col)
    return lines



def calcPos(f, e, i, a, w, om):
    i = math.radians(i)
    w = math.radians(w)
    om = math.radians(om)
    p = a* (1 - e**2)
    tan = np.tan
    sini = sin(i)
    cosi = cos(i)
    r = p/(1 + e *cos(f))
    sinb = sini * math.sin(w + f)
    tanlom = cosi * tan(w + f)
    lom = np.arctan(tanlom)
    if (sin(lom)* sin(w+f)) <  0:
        lom += pi
    l = lom + om
    b = np.arcsin(sinb)

    x = r * cos(b)*sin(l)
    y = r * cos(b)*cos(l)
    z = r * sin(b)
    return [x, y, z]

ax.scatter([0], [0], [0], color = 'yellow')
ax.scatter([0], [0], [-5], color = 'white', alpha = 0)
ax.scatter([0], [0], [5], color = 'white', alpha = 0)
ax.tick_params(axis ='x', length = 1.0)
ax.tick_params(axis ='z', length = 1.0)
ax.tick_params(axis ='y', length = 1.0)

ax.grid(False)
plt.axis('off')
#for vesta:
#x1, y1, z1 = calcPos(f, e, i, a, w, om)
#ax.scatter([x1], [y1], [z1], color = 'r')
#for earth:
#x2, y2, z2 = calcPos(0, 0.0167086, 0.00005, 1, 114.21, -11.26064)
#ax.scatter([x2], [y2], [z2], color = 'blue')
#
#for mars:
#x3, y3, z3 = calcPos(0, 0.0934, 1.850, 1.5237, 286.5, 49.57854)
#ax.scatter([x3], [y3], [z3], color = 'orange')
#
par1 = [0, 0.0167086, 0.00005, 1, 114.21, -11.26064]
par2 = [0, e, i, a, w, om]
par3 = [0, 0.0934, 1.850, 1.5237, 286.5, 49.57854] 

#img = plt.imread("backgr.jpeg")
inc0 = 2 * pi /1000
#ax.imshow(img)
prepareOrbit(f, e, i, a, w, om, 'g')
prepareOrbit(0, 0.0167086, 0.00005, 1, 114.21, -11.26064, 'b')
prepareOrbit(0, 0.0934, 1.850, 1.5237, 286.5, 49.57854, 'r')
prepareOrbit(0, 0.0489, 1.303, 5.2038, 273.867, 100.464, 'k')

x1, y1, z1 = precalcPos(e, i, a, w, om, 10, 0)
x2, y2, z2 = precalcPos(0.0167086, 0.00005, 1, 114.21, -11.26064, 10, 0)
x3, y3, z3 = precalcPos(0.0934, 1.850, 1.5237, 286.5, 49.57854, 10, 0)
x4, y4, z4 = precalcPos(0.0489, 1.303,5.2038, 273.867, 100.464, 10, 0)
#print(x)
#print(y)
#print(z)
ax.set_title('Visualisation of orbits')
#ax.scatter3D(x, y, z, cmap='cividis')
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)
#print(x1)
#print(y1)
#print(z1)

myAnimation = animation.FuncAnimation(fig, animate, frames = 1500, interval = 1, repeat=True)
#myAnimation.save('orbits.gif')
plt.show()


    



