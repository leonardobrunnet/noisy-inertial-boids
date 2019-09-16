#Details:
#Particles have different sizes (Req=1+a*rand()) a~0.1
#Particles reaching the right edge are reinserted in the first left quarter with the same
#speed v, direction n, and y-position.
#good noise value 1.0
#good density 1.0
#good value for wall_osc~0.01

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import math as math
import random as rand

#Particle class definition
class particle:
   # noise=0.6 #original value
    noise_par=1.0
    noise_per=0.1
    gamma=0.1
    rq = 0.5
    v0 = 1.0
    mu=1.0
    Frep=30.0
    Fadh=0.75
    #Req=5./6. # original work by Szabo
    #R0=1.0
#    Req=1.0
    R0=6./5.
    def __init__(self, x, y, vx, vy, ident, Raio_equilibrio):
        self.r = np.array([x,y])
        self.v =  np.array([vx*self.v0,vy*self.v0])
        self.theta = np.arctan2(vy,vx)
        self.n_par = np.array([np.cos(self.theta),np.sin(self.theta)])
        self.n_per = np.array([np.sin(self.theta),-np.cos(self.theta)])
        self.v_par = np.dot(self.v,self.n_par)
        self.ident = ident
        self.Mybox = int((self.r[0]+L[0])/lbox)+nb[0]*int((self.r[1]+L[1])/lbox)
        self.Force =np.array([0.,0.])
        self.Req = 1+Raio_equilibrio
        
    
    def mybox(self): #Each particle calculates the box it is in
        j=int((self.r[0]+L[0])/lbox)+nb[0]*int((self.r[1]+L[1])/lbox)
        return j

    def changebox(self):
        newbox=self.mybox()
        if(newbox!=self.Mybox): #verify particle box change 
            box[self.Mybox].mylist.remove(self.ident) #this is safe since particles ident is unique
            box[newbox].mylist.append(self.ident)
            self.Mybox=newbox
        return self.Mybox

    def mov(self): #Particle moviment
        self.v_par+=-self.gamma*self.v_par*dt+self.noise_par*(rand.random()-0.5)*np.sqrt(dt)
        dr_parallel=self.v_par*dt
        dr_perpend = noise_per*(rand.random()-0.5)*np.sqrt(dt)
        dr=dr_par*self.n+dr_per*self.n_per
        self.r+=dr
        self.theta+=self.rq*noise_per*(rand.random()-0.5)*sqrt(dt)
        self.n_par = np.array([np.cos(self.theta),np.sin(self.theta)])
        self.n_per = np.array([np.sin(self.theta),-np.cos(self.theta)])
        self.contour()
        return self.r,self.v, self.n

    
    def contour(self):
        if self.r[0]<-L[0]:
            self.n[0]=np.abs(self.n[0])
            self.v[0]=np.abs(self.v[0])
            self.theta=(rand.random()-0.5)*3.14
#            self.r[0]=-2*L[0]-self.r[0]
            self.r[0]=-L[0]+wall_osc*rand.random()*self.Req
        if self.r[0]>L[0]: #will be send back to the first quarter, v, n and theta remain the same
            self.r[0]=(-1+rand.random()/4.)*L[0]
        if self.r[1]<-L[1]:
            self.n[1]=np.abs(self.n[1])
            self.v[1]=np.abs(self.v[1])
            self.theta=rand.random()*3.14
            self.r[1]=-L[1]+wall_osc*self.Req*rand.random()
        if self.r[1]>L[1]:
            self.n[1]=-np.abs(self.n[1])
            self.v[1]=-np.abs(self.v[1])
            self.theta=-rand.random()*3.14
            self.r[1]=L[1]-wall_osc*self.Req*rand.random()
        normr=np.linalg.norm(self.r)
        #Particles hitting the cylinder loose their velocity component perpendicular to the cylinder if invading it
        # if normr < cylinder_radius:
        #     vaux=self.v
        #     vpar=np.dot(self.v,self.r)*self.r/np.linalg.norm(self.r)**2 #velocity component parallel to cylinder radius
        #     vper=self.v-vpar #perpendicular component to cylinder radius (parallel to cylinder contour)
        #     crit = np.dot(vpar,self.r)
        #     if crit < 0 :
        #         self.v=vper
        #         self.r-=vaux*dt
        #     self.r+=self.v*dt 
        #     self.n=self.v/np.linalg.norm(self.v)
        #     self.theta = np.arctan2(self.v[1],self.v[0])
        return self.r, self.theta, self.n, self.v

    def autovelchange(self,v,n):
        vnorm=np.linalg.norm(v)
        u=self.v/vnorm
        self.dtheta=np.arcsin(self.n[0]*u[1]-self.n[1]*u[0])
        return self.dtheta

    def forces_between_particles(self):
        def force(self,dr,Req):

            normdr=np.linalg.norm(dr)
            
            if(normdr<Req):
                f=self.Frep*(1/self.Req-1/normdr)*dr
            else:
                if(normdr<self.R0):
                    f=self.Fadh*dr*(1-Req/normdr)/(self.R0-Req)
                else:
                    f=0.0
            return f
        for i in box[self.Mybox].mylist:
            if(self.ident!=i):
                dr=part[i].r-self.r
                Req=(part[i].Req+self.Req)/2.
                #print self.ident, i, np.linalg.norm(dr)
                self.Force+=force(self,dr,Req)
        for i in box[self.Mybox].neighboxlist:
            Req=(part[i].Req+self.Req)/2.
            dr=part[i].r-self.r
            #            print self.ident, i, np.linalg.norm(dr),t
            f=force(self,dr,Req)
            self.Force+=f
            part[i].Force-=f
    def zeroforce(self):
        self.Force=np.array([0.,0.])
        

#Box class definition
class boite:
    def __init__(self,index):
        self.index = index
        self.mylist = []
        self.neighboxlist = []

    def neighbor_list(self):
        self.neighboxlist=[]
        if self.index%nb[0]!=nb[0]-1:
            self.neighboxlist.extend(box[self.index+1].mylist)
        if self.index%nb[0]!=0 and self.index+nb[0]<nb2:
            self.neighboxlist.extend(box[self.index+nb[0]-1].mylist)
        if self.index+nb[0]<nb2:
            self.neighboxlist.extend(box[self.index+nb[0]].mylist)
        if self.index+nb[0]<nb2 and self.index%nb[0]!=nb[0]-1:
            self.neighboxlist.extend(box[self.index+nb[0]+1].mylist)

#Main program                                  
#global variables
global N,L,lbox,nb,dt,nb2,t,cylinder_radius
N=500
L=np.array([15,5])
lbox=1
nb=2*L
nb2=nb[1]*nb[0]
dt=0.01
exit_fig=20
tau=10.0
passos=50000
rand.seed(0.1)
cylinder_radius=2.25
wall_osc = 0.01  #amplitude of random number separating a boid from the walls
size_disp = 0.1  #particle size Req dispersion
output_file_name="output_simu_sabo.dat"
output_file = open(output_file_name,'w')
output_file.write("Number_of_particles: %d\n"%N)
output_file.write("Box_size: %d\n"%lbox)
output_file.write("Steps_between_images: %d \n"%exit_fig)
output_file.write("Radius: %f \n"%cylinder_radius)
output_file.write("Dimensions: %d %d \n"%(L[0],L[1]))
output_file.write("Max-dist: 4 \n")
                  
#initialize N particles

#part=list(particle(-L[0]+L[0]/2*rand.random(),L[1]*2*(rand.random()-0.5), rand.random()-0.5,rand.random()-0.5, i) for i in range(N))
part=list(particle(-L[0]+3*L[0]/4.*rand.random(),L[1]*2*(rand.random()-0.5), rand.random()-0.5,rand.random()-0.5, i, size_disp*(rand.random()-0.5) ) for i in range(N))
box=list(boite(i) for i in range(nb2))
t=0

# Construct list of particles in each box

for i in range(N):
    part[i].Mybox=part[i].mybox()
    box[part[i].Mybox].mylist.append(i)

#print part[98].Mybox,box[part[98].Mybox].mylist


    # Construct list of particles in neighboring boxes

map(lambda i:i.neighbor_list(), box)


#print box[part[98].Mybox].neighboxlist


#System evolution
figindex=0
intt=0
while(t<passos*dt):

#Calculate the forces
    map(lambda i:i.forces_between_particles(), part)
#Move all particles
    map(lambda i:i.mov(), part)
    t+=dt #update time
#Find the newboxes 
    map(lambda i:i.changebox(), part)
#Construct the list of particles in neighboring boxes
    map(lambda i:i.neighbor_list(), box)
#Reset forces to zero
    map(lambda i:i.zeroforce(), part)

# Make a scatter graph
    intt+=1
    delta=5.
    sizes=1.5
    if(intt%exit_fig==0):
        print(t)
        output_file.write("x y vx vy \n")
        plt.axis([-L[0]-delta,L[0]+delta,-L[1]-delta,L[1]+delta])
        plt.axes().set_aspect(1.0)
        circle=plt.Circle((0.,0.),radius=cylinder_radius-0.5,color='r')
        x,y,vx,vy=[],[],[],[]
        map(lambda i:x.append(i.r[0]), part)
        map(lambda i:y.append(i.r[1]), part)
        map(lambda i:vx.append(i.v[0]), part)
        map(lambda i:vy.append(i.v[1]), part)
        for i in range(len(x)):
            output_file.write("%f %f %f %f \n"%(x[i],y[i],vx[i],vy[i]))

        plt.scatter(x,y,s=sizes,alpha=0.3)
        name=str(figindex)+".png"
        fig = plt.gcf()
        plt.rc("savefig",dpi=300)
        ax = fig.gca()
#        ax.add_artist(circle)
        fig.savefig(name,bbox_inches='tight')
        figindex+=1
        fig.clf()
