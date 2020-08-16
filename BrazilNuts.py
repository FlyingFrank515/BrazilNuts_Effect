from vpython import *
import numpy as np
N = 200
m, size = 1, 0.5 #1 gram and the radius is 0.5cm
L = 15 # 2L is the cubic container's original length, width, and height
t, dt = 0, 0.001
g=980 #cm/s^2
nuts = [] # list to store atoms
scene = canvas(width=500, height=500, background=vector(0.2,0.2,0), align = 'left')
container = box(length = 2*L, height = 2*L, width = 2*L, opacity=0.2, color = color.yellow )
p_a, v_a = np.zeros((N,3)), np.zeros((N,3)) # particle position array and particle velocity array, N particles and 3 for x, y, z
for i in range(N):
    p_a[i] = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L] # particle is initially random positioned in container
    nut = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=color.orange)
    ra = pi*random()
    rb = 2*pi*random()
    v_a[i] = [0,0,0]
    nuts.append(nut)
p_bignut=np.zeros((N,3))
p_bignut = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L]
bignut = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = 2, color=color.red)
v_bignut=[0,0,0]
def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime
while True:
    t += dt
    rate(1000)
    p_a += v_a*dt-g*dt # calculate new positions for all atoms
    p_bignut+=v_bignut*dt-g*dt
    for i in range(N): nuts[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
    bignut.pos=vector(p_bignut[i,0],p_bignut[i,1],p_bignut[i,2])
    ### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis] # array for vector from one atom to another atom for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) # distance array between atoms for all pairs of atoms
    hit = np.less_equal(rmag,2*size)-np.identity(N) # if smaller than 2*size meaning these two atoms might hit each other
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() # change hit to a list
    for ij in hitlist: # i,j encoded as i*Natoms+j
        i, j = divmod(ij,N) # atom pair, i-th and j-th atoms, hit each other
        hitlist.remove(j*N+i) # remove j,i pair from list to avoid handling the collision twice
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 : # only handling collision if two atoms are approaching each other
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision
    #find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        if abs(p_a[i][0]) >= L - size and p_a[i][0]*v_a[i][0] > 0 :
            v_a[i][0] = - v_a[i][0]
        if abs(p_a[i][1]) >= L - size and p_a[i][1]*v_a[i][1] > 0 :
            v_a[i][1] = - v_a[i][1]
        if abs(p_a[i][2]) >= L - size and p_a[i][2]*v_a[i][2] > 0 :
            v_a[i][2] = - v_a[i][2]