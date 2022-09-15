import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.animation as animation
from scipy.integrate import odeint

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
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
gm = 39.4784176
def drdt(u, t):
    temp = -gm/((u[0]**2+u[1]**2+u[2]**2)**1.5)
    return (u[3], u[4], u[5], temp*u[0], temp*u[1], temp*u[2])
   

def angle(sinx, cosx):
    return math.atan2(sinx, cosx)

def calcE(M, e):
    M = M%(2 * pi)
    E = M
    for i in range(8):
        E = M + e * sin(E) 
    return E

def updateRF(a, e, t, Per):
    Per = Per * 365.25
    M = 6.2831853/Per * t
    E_eSinE = M
    E = calcE(M, e)
    cosE = cos(E)
    rnow = a*(1 - e *cosE)
    f = 2*np.arctan(math.sqrt((1+e)/(1-e)*tan(E/2)))
    return (rnow, f)

def drawOrbit(e, i, a, w, om, col):
    e, i, a, w, om = map(float, [e, i, a, w, om])
    i, w, om = map(math.radians, [i, w, om])
    p = a* (1 - e**2)
    sini = sin(i); cosi = cos(i)
    f = 0
    inc = 2 *pi/1000
    x, y, z = [np.empty(1000) for i in range(3)]
    for i in range(1000):
        cosf = cos(f)
        r = p/(1 + e*cosf)
        sinb = sini * sin(w + f)
        tanlom = cosi * tan(w + f)
        lom = np.arctan(tanlom)
        if (sin(lom)* sin(w+f)) <  0:
            lom += pi
        l = lom + om
        b = arcsin(sinb)

        x[i] = r * cos(b)*sin(l)
        y[i] = r * cos(b)*cos(l)
        z[i] = r * sin(b)
        f += inc

    ax.plot(x, y, z, col)
    
def updateVect(pos, vel, dt):
    temp = gm/((pos[0]**2 + pos[1]**2+pos[2]**2)**1.5)
    acc = -temp * pos
    pos += vel * dt
    vel += acc * dt
    return (pos, vel)

def initPos(r, a, i, w, om):
    cosi = cos(i); sini = sin(i)
    pos = np.empty(3)
    sinw = sin(w); cosw= cos(w)
    sinom = sin(om); cosom=cos(om)
    sinb = sini * math.sin(w)
    tanlom = cosi * tan(w)
    lom = np.arctan(tanlom)
    if (sin(lom)* sin(w)) <  0:
        lom += pi
    l = lom + om
    b = np.arcsin(sinb)

    pos[0] = r * cos(b)*sin(l)
    pos[1] = r * cos(b)*cos(l)
    pos[2] = r * sin(b)
    vel = np.empty(3)
    velMag = math.sqrt(gm*(2/r - 1/a))
    vel[1] = (-cosom*sinw - sinom*cosw*cosi)*velMag
    vel[0] = (-sinom*sinw+cosom*cosw*cosi)*velMag
    vel[2] = cosw*sini*velMag
    return pos, vel
        
def calcPos(e, i, a, w, om, dt):
    e, i, a, w, om, dt = map(float, [e, i, a, w, om, dt])
    i, w, om = map(math.radians, [i, w, om])
    Per = a**1.5
    sini = sin(i); cosi = cos(i)
    t = 0
    dt = dt/365.25
    r, f = updateRF(a, e, t, Per)
    pos, vel = [np.empty(0) for _ in range(2)]
    pos, vel = initPos(r, a, i, w, om)    
    
    tarr = np.arange(0, 10, 10/animFrames)
    y0 =  [*pos, *vel]
    loc = odeint(drdt, y0, tarr)
    x = loc[:,0]; y = loc[:,1];z=loc[:,2]
    
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

#ax.grid(False)
#plt.axis('off')	

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
