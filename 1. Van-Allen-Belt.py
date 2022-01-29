from visual import*
from random import random

N = 500
m,size = 1.7E-27,100000
SIZE = 6371E3             #Earth radius
B0 = 3.12E-5
vrms = 450E4
atoms = []
q = 1.6E-19         #bring charge
R_E = 6372000       #radius of earth
t,dt = 0,0.001

scene = display(width=800,height=800,center=(-R_E*2,0,0),background=(0,0.2,0))

def MFa(v,pos_array):    #magnetic force's accelerlator
    theta = arccos(dot(pos_array, (0,1,0))/mag(pos_array))
    r = mag(pos_array)
    c = norm(cross(cross((0, 1, 0), pos_array), pos_array))
    B = -2*B0*((R_E/r)**3)*cos(theta)*norm(pos_array) - B0 * ((R_E/r) ** 3) * sin(theta) * norm(c) 
    F=q*cross(v,B)
    a=F/m
    return a

pos_array,v_array,a_array = zeros((N,3)),zeros((N,3)),zeros((N,3))
earth = sphere(pos = (0,0,0),radius=SIZE,material=materials.earth)
for i in range(N):
    ra = 2*pi*random()
    rb = 2*pi*random()
    pos_array[i] = [-(R_E)*2,2*R_E*cos(ra)*sin(rb),2*R_E*sin(ra)*sin(rb)]
    atom = sphere(pos = pos_array[i],radius=size,color=color.yellow,make_trail=True,retain=10)
    v_array[i] = [vrms*100,0,0]
    atoms.append(atom)

while True:
    t+=dt
    rate(1000)
    for i in range(N):
        a_array[i] = MFa(v_array[i],pos_array[i])
    atoms.append(atom)

    pos_array+= v_array*dt
    v_array+= a_array*dt
    for i in range(N):
        atoms[i].pos = pos_array[i]
        atoms[i].v = v_array[i]
        
        if mag(pos_array[i])>= R_E*10:
            ra = 2*pi*random()
            rb = 2*pi*random()
            pos_array[i] = [-(R_E)*2,2*R_E*cos(ra)*sin(rb),2*R_E*sin(ra)*sin(rb)]
            v_array[i] = [vrms*100,0,0]