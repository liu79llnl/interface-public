import math

def translation(vx, vy, dt):
    velocity = lambda p : [dt*vx, dt*vy]
    return velocity

#Rotation counterclockwise about center at 1 radian per second
def rotation(center, dt):
    velocity = lambda p : [-dt*(p[1]-center[1]), dt*(p[0]-center[0])]
    return velocity

#Vortex example from Rider and Kothe
def vortex(dt, time_reverse):
    velocity = lambda p : [-2*dt*(math.sin(math.pi*p[0]))**2 * math.sin(math.pi*p[1]) * math.cos(math.pi*p[1]), 2*dt*math.sin(math.pi*p[0])*math.cos(math.pi*p[0]) * (math.sin(math.pi*p[1]))**2]
    return velocity
    
#Reverses at time time_reverse (after time_reverse/dt time steps)
def vortex_tr(dt, t, totalt):
    velocity = lambda p : [-2*dt*math.cos(math.pi*t/totalt)*(math.sin(math.pi*p[0]))**2 * math.sin(math.pi*p[1]) * math.cos(math.pi*p[1]), 2*dt*math.cos(math.pi*t/totalt)*math.sin(math.pi*p[0])*math.cos(math.pi*p[0]) * (math.sin(math.pi*p[1]))**2]
    return velocity