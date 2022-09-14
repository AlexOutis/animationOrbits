import numpy as np
import matplotlib.pyplot as plt
import math 
from mpl_toolkits import mplot3d
import matplotlib.animation as animation

fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(111, projection='3d')
cos = np.cos;sin = np.sin;tan = np.tan
arcsin = np.arcsin;arccos = np.arccos; pi = math.pi; floor = math.floor
params = np.empty((0, 6))
filename = input("Input file's name: ")
f = open(filename, 'r')

def addParams(params, ln):
    e, i, a, w, om = map(float, ln.split()[:5])
    col = ln.split()[5]
    params = np.append(params, np.array([[e, i, a, w, om, col]]), axis=0)
    return params
ln = f.readline()

while ln:
    params = addParams(params, ln)
    ln = f.readline()


dots = []

years = int(input("Choose animation duration (in years, for example 10): "))
animDur = math.floor(years * 365.25)
animSpeed = 10
animFrames = math.floor(animDur / animSpeed)
def angle(sinx, cosx):
    if sinx > 0 and  cosx > 0:
        return arcsin(sinx)
    elif sinx > 0 and cosx < 0:
        return arccos(cosx)
    elif sinx < 0 and cosx < 0:
        return pi -  arcsin(sinx)
    else: return arcsin(sinx)

def calcE(M, e):
    E = M
    M = M%(2 * pi)
    for i in range(10):
        E = M + e * sin(E) 
    return E

def updateRF(a, e, t, Per):
    Per = Per * 365.25
    M = 2 * pi/Per * t
    E_eSinE = M
    E = calcE(M, e)
    rnow = a*(1 - e *cos(E))
    cosf = (cos(E) - e)/(1 - e * cos(E))
    sinf = math.sqrt(1 - e**2) * sin(E)/(1 - e * cos(E))
    f = angle(sinf, cosf) 
    return (rnow, f)

def drawOrbit(e, i, a, w, om, col):
    e = float(e); i = float(i); a= float(a); w = float(w); om = float(om); 
    i = math.radians(i)
    w = math.radians(w)
    om = math.radians(om)
    p = a* (1 - e**2)
    Per = math.sqrt(a**3)
    sini = sin(i)
    cosi = cos(i)
    f = 0
    inc = 2 *pi/1000
    x, y, z = [np.empty(1000) for i in range(3)]
    for i in range(1000):
        cosf = cos(f)
        r = p/(1 + e*cosf)
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

    ax.plot(x, y, z, col)
    



def calcPos(e, i, a, w, om, dt):
    e = float(e); i = float(i); a= float(a); w = float(w); om = float(om); dt = float(dt)
    i = math.radians(i)
    w = math.radians(w)
    om = math.radians(om)
    p = a* (1 - e**2)
    Per = math.sqrt(a**3)
    sini = sin(i)
    cosi = cos(i)
    t = 0
    x, y, z = [np.empty(animFrames) for i in range(3)]
    for i in range(animFrames):
        r, f = updateRF(a,e, t, Per)
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
    return [x, y, z]
    
def init_dots():
    for param in params:
        dots.append(ax.plot([0], [0], [0], f'{param[5]}o')[0])

def animate(k):
    ax.set_title(f'{k*animSpeed} days')
    for j in range(params.shape[0]):
        dots[j].set_data_3d(coords[j][0][k], coords[j][1][k], coords[j][2][k])
    
    return dots

def maxLim():
    x = 0
    for param in params:
        x = max(float(param[2]), x)
    return x

ax.plot([0], [0], [0], 'yo')
maxDist = maxLim()*1.3
ax.set_xlim(-maxDist, maxDist)
ax.set_ylim(-maxDist, maxDist)
ax.set_zlim(-maxDist, maxDist)

ax.grid(False)
plt.axis('off')

init_dots()

coords = np.empty([0, 3, animFrames])

for param in params:
    drawOrbit(*param)
    temp = np.array(calcPos(*param[:5], animSpeed))
    temp = np.reshape(temp, (1, 3, animFrames))
    coords = np.append(coords, temp, axis=0)

ax.set_title('Visualisation of orbits')

myAnimation = animation.FuncAnimation(fig, animate, frames = animFrames, interval = 40, repeat=True)
plt.show()
