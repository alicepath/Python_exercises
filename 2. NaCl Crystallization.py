from visual import*
from random import random
N = 50
k, T = 1.38E-23, 50
t, dt = 0, 0.5E-13
q = 1.6 * 10 ** -19
u = 1.66 * 10 ** -27
mn = 22.9897693 * u                     #mass of Na
mc = 35.45 * u                          #mass of Cl
L = ((24.4E-3/(6E23)) * N) ** (1/3.0)/7 #size of container
sizeN = 227E-12                         #size of Na
sizeCl = 175E-12                        #size of Cl
L_sizeN = L - sizeN
L_sizeCl = L - sizeCl
vrms_Na = 1E-6*(3*k*T/mn)**0.5
vrms_Cl = 1E-6*(3*k*T/mc)**0.5
ep = 1E5 * 8.854E-12

atomN = []
atomCl = []
scene = display(width = 800, height = 800, background = (0.2, 0.2, 0))
container = box(length = 2 * L, height = 2 * L, width = 2 * L, opacity = 0.2, color = color.yellow)


def v_collision(a1p,a2p,a1v,a2v,a1m,a2m):              #function after collision velocity
    v1_prime=0.9 * (a1v-2*a2m/(a1m+a2m)*(a1p- a2p) *sum((a1v-a2v)*(a1p-a2p))/sum((a1p-a2p)**2))
    v2_prime=0.9 * (a2v-2*a1m/(a2m+a1m)*(a2p- a1p) *sum((a2v-a1v)*(a2p-a1p))/sum((a2p-a1p)**2))
    return(v1_prime,v2_prime)

Npos_array,Nv_array, Na_array1,Na_array2,Na_array = zeros((N,3)),zeros((N,3)), zeros((N,3)), zeros((N,3)), zeros((N,3))
for i in range(N):

    Npos_array[i] = [-L_sizeN + 2*L_sizeN*random(), -L_sizeN + 2*L_sizeN*random(), -L_sizeN + 2*L_sizeN*random()]
    ra, rb = pi*random(), 2*pi*random()
    Nv_array[i] = [vrms_Na*sin(ra)*cos(rb), vrms_Na*sin(ra)*sin(rb), vrms_Na*cos(ra)]
    atom = sphere(pos = Npos_array[i],radius=sizeN,color=color.white)
    atomN.append(atom)

Clpos_array,Clv_array, Cla_array1, Cla_array2, Cla_array = zeros((N,3)),zeros((N,3)), zeros((N,3)), zeros((N,3)), zeros((N,3))
for i in range(N):
    Clpos_array[i] = [-L_sizeCl + 2*L_sizeCl*random(), -L_sizeCl + 2*L_sizeCl*random(), -L_sizeCl + 2*L_sizeCl*random()]
    ra, rb = pi*random(), 2*pi*random()
    Clv_array[i] = [vrms_Cl*sin(ra)*cos(rb), vrms_Cl*sin(ra)*sin(rb), vrms_Cl*cos(ra)]
    atom = sphere(pos = Clpos_array[i],radius = sizeCl, color = color.green)    
    atomCl.append(atom)
    
    

while 1:
    rate(100000)

    r1_array = Npos_array-Npos_array[:,newaxis]                #all pairs of atom-to-tom vectors
    r2_array = Clpos_array-Clpos_array[:,newaxis]
    r3_array = Npos_array-Clpos_array[:,newaxis]

    r1mag = sqrt(sum(square(r1_array),-1))                    #atom-to-atom scalar distances
    r2mag = sqrt(sum(square(r2_array),-1))
    r3mag = sqrt(sum(square(r3_array),-1))
    QQ1 = greater(r1mag,0.000000000000000000000001)
    QQ2 = greater(r2mag,0.000000000000000000000001)
    QQ3 = greater(r3mag,sizeN + sizeCl)

    hitlistQ1 = sort(nonzero(QQ1.flat)[0]).tolist()
    hitlistQ2 = sort(nonzero(QQ2.flat)[0]).tolist()
    hitlistQ3 = sort(nonzero(QQ3.flat)[0]).tolist()
    for ij in hitlistQ1 :
        i,j = divmod(ij,N)
#        hitlistQ1.remove(j*N+i)
        Na_array1[i] +=  -(1/(4 * pi * ep) * q ** 2 / (r1mag[i][j])**3 * (r1_array[i][j]))/mn
    for ij in hitlistQ2 :
        i,j = divmod(ij,N)
#        hitlistQ1.remove(j*N+i)
        Cla_array1[i] += -(1/(4 * pi * ep) * q ** 2 / (r2mag[i][j])**3 * (r2_array[i][j]))/mc

    for ij in hitlistQ3 :
        i,j = divmod(ij,N)
#        hitlistQ1.remove(j*N+i)
        Na_array2[i] += -(1/(4 * pi * ep) * q ** 2/(r3mag[i][j])**3 * (r3_array[i][j]))/mn
        Cla_array2[i] += (1/(4 * pi * ep) * q ** 2/(r3mag[i][j])**3 * (r3_array[i][j]))/mc
    Na_array = (Na_array1 + Na_array2)
    Cla_array = (Cla_array1 + Cla_array2)

    
    Nv_array += Na_array * dt    
    Npos_array+= Nv_array*dt
    Clv_array += Cla_array * dt
    Clpos_array+= Clv_array*dt
    for i in range(N):
        atomN[i].pos = Npos_array[i]
        atomCl[i].pos = Clpos_array[i]
#######################################################################################################################################
#find collisions between pairs of atoms, and handle their collisions
    hit1 = less_equal(r1mag,2*sizeN)-identity(N)               #find out those atom-to-atom distances smaller than 2*size
    hit2 = less_equal(r2mag,2*sizeCl)-identity(N)
    hit3 = less_equal(r3mag,sizeN+sizeCl)

    hitlist1 = sort(nonzero(hit1.flat)[0]).tolist()           #i,j encoded as i*Natoms+j
    for ij in hitlist1 :
        i,j = divmod(ij,N)                                  #decode atom pair
        hitlist1.remove(j*N+i)                               #remove symmetric j,i pair from list
        if sum((Npos_array[i]-Npos_array[j])*(Nv_array[i]-Nv_array[j]))<0:
            Nv_array[i],Nv_array[j] = v_collision(Npos_array[i],Npos_array[j],Nv_array[i],Nv_array[j],mn,mn)
    
    hitlist2 = sort(nonzero(hit2.flat)[0]).tolist()           #i,j encoded as i*Natoms+j
    for ij in hitlist2 :
        i,j = divmod(ij,N)                                  #decode atom pair
        hitlist2.remove(j*N+i)                               #remove symmetric j,i pair from list
        if sum((Clpos_array[i]-Clpos_array[j])*(Clv_array[i]-Clv_array[j]))<0:
            Clv_array[i],Clv_array[j] = v_collision(Clpos_array[i],Clpos_array[j],Clv_array[i],Clv_array[j],mc,mc)

    hitlist3 = sort(nonzero(hit3.flat)[0]).tolist()           #i,j encoded as i*Natoms+j
    for ij in hitlist3 :
        i,j=divmod(ij,N)                 #decode atom pair                          
        if sum((Npos_array[i]-Clpos_array[j])*(Nv_array[i]-Clv_array[j]))<0:
            Nv_array[i],Clv_array[j] = v_collision(Npos_array[i],Clpos_array[j],Nv_array[i],Clv_array[j],mn,mc)
#find collisions between the atoms and the walls, and handle their collisions
    for i in range(N):
        if abs(Npos_array[i][0]) >= L_sizeN and Npos_array[i][0]*Nv_array[i][0] > 0:
            Nv_array[i][0]= - Nv_array[i][0] - Na_array[i][0] * dt
        if abs(Npos_array[i][1]) >=L_sizeN and Npos_array[i][1]*Nv_array[i][1] > 0:
            Nv_array[i][1]= - Nv_array[i][1] - Na_array[i][1] * dt
        if abs(Npos_array[i][2]) >=L_sizeN and Npos_array[i][2]*Nv_array[i][2] > 0:
            Nv_array[i][2]= - Nv_array[i][2] - Na_array[i][2] * dt

        if abs(Clpos_array[i][0]) >= L_sizeCl and Clpos_array[i][0]*Clv_array[i][0] > 0:
            Clv_array[i][0]=- Clv_array[i][0] - Cla_array[i][0] * dt
        if abs(Clpos_array[i][1]) >= L_sizeCl and Clpos_array[i][1]*Clv_array[i][1] > 0:
            Clv_array[i][1]=- Clv_array[i][1] - Cla_array[i][1] * dt
        if abs(Clpos_array[i][2]) >= L_sizeCl and Clpos_array[i][2]*Clv_array[i][2] > 0:
            Clv_array[i][2]=- Clv_array[i][2] - Cla_array[i][2] * dt