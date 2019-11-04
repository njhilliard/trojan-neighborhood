#!/usr/bin/env python
# coding: utf-8

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import readsnap as rs
#import readsnapHDF5 as rs
import h5py 
import pylab
import struct
#import time
#import readhsml
import matplotlib.cm as cm
import conversions as co
#from PIL import Image
#import PIL.ImageOps    

#parttype=0
parttype=2
arepo=0


#base='/home/donghia/MGroups/output_arms/'
#base='/home/donghia/MWbarnogas/Run/medres/output/'
base='./data/'

tname=" "
#output_dir="/home/donghia/MWbarnogas/Run/medres/output/movie/faceon/"
output_dir="./analysis/"
#base='/n/hernquistfs1/mvogelsberger/projects/spirals/perturbers/RUN_c9_mc_pro_finite_short/output/'
#tname="MCs prograde short"
#output_dir="./output/c9_mc_pro_finite_short/"

#base='/n/hernquistfs1/mvogelsberger/projects/spirals/perturbers/RUN_c9_mc_pro_finite_long/output/'
#tname="MCs prograde long"
#output_dir="./output/c9_mc_pro_finite_long/"

#base='/n/hernquistfs1/mvogelsberger/projects/spirals/perturbers/RUN_c9_mc_pro_turnMCsoff/output/'
#tname="MCs turned off"
#output_dir="/n/hernquistfs1/edonghia/Spirals/Movie_turnoff/"
#output_dir="/n/scratch2/hernquist_lab/edonghia/Movie_swithoff_M101/"
#output_dir="/n/scratch2/hernquist_lab/edonghia/SwichallGMCtime20/output/MovieSwitch20/"
#output_dir="/n/itcbackup1/edonghia/MovieSwichallGMCtime10/"

snapbase='snap_'

#list of snapshots
nums=[582]


#plot does for MCs (1=on)
mc_flag=0
#size of MC dots
mc_s=2.5
patch_flag=1
patch_s=0.5



#select axes
xax=0
yax=1
zax=2

#for mass weighted histogram
BINS=512

#for deviation field
BINS_theta=360
BINS_r=360


#for fitting exponential model
rmin=0.0
rmax=15.0
#Rs=2.8

#bin points
dtheta=2*np.pi/BINS_theta
dr=(rmax-rmin)/BINS_r
thetamid=(np.arange(BINS_theta)+0.5) * dtheta - np.pi   
rmid=(np.arange(BINS_r)+0.5) * dr + rmin


#size in units of scale length
lengthX=15.0
lengthY=15.0
vx0=-5.0
vy0=170.0
#lengthX=300.0-vx0
#lengthY=300.0-vy0


Zmin=-4.25
Zmax=-0.27
#Zmin=0.0
#Zmax=0.0

urnew=np.zeros(15000000)
utnew=np.zeros(15000000)
pznew=np.zeros(15000000)

px=np.zeros(15000000)
py=np.zeros(15000000)
pz=np.zeros(15000000)

for i in range(0,len(nums)):
        num=nums[i]

        #####commented: how to read in hdf5
        filename=base+snapbase+str(num).zfill(3)+".hdf5"        
        #filename=base+"snapdir_"+str(num).zfill(3)+"/"+snapbase+str(num).zfill(3)
        print(filename)
        with h5py.File(filename,"r") as filename:
        #with h5py.File("/home/donghia/MW_model_x20_f20_VC16hd_new_snapshot_306.hdf5","r") as filename:

                ###this indented block must be done to access all data inside the hdf5 file opened as 'filename'
                myHead = filename.__getitem__("/Header")
                pDisk = filename.__getitem__("/PartType2")
                pBulge = filename.__getitem__("/PartType3")
                diskCoord = np.array(pDisk.__getitem__("Coordinates"))
                diskVels = np.array(pDisk.__getitem__("Velocities"))
                bulgeCoord = np.array(pBulge.__getitem__("Coordinates"))
                bulgeVels = np.array(pBulge.__getitem__("Velocities"))
                ####your coords and vels are arrays of size (nParticles_disk,3)
                att = myHead.attrs  ###your header attributes
                myTime = att.__getitem__("Time")
                myMasses = np.array(att.__getitem__("MassTable"))
                myDiskMass = myMasses[2]
                myBulgeMass = myMasses[3] 
                print( 'Processing ', filename          )

        ###here after the indentation of the 'with' block is done, python will close the hdf5 file correctly and safe

        #head=rs.snapshot_header(filename)
        #print( "time = ", head.time)
        print( "time = ", myTime)
        print( "length of coordinate of disk is ", diskCoord.shape)
        print( "length of coordinate of bulge is ", bulgeCoord.shape)
        print( "myDiskMass is ", myMasses[2], myMasses[3])


        #pos=rs.read_block(filename, "POS ", parttype=2).astype('float64')
        #poss=rs.read_block(filename, "POS ", parttype=2).astype('float64')
        #vel=rs.read_block(filename, "VEL ", parttype=2).astype('float64')
        ##pos=rs.read_block(filename, "POS ", parttype=3).astype('float64')
        #mass=rs.read_block(filename, "MASS", parttype=2).astype('float64')
        ##mass=rs.read_block(filename, "MASS", parttype=3).astype('float64')
           

        #px=pos[:,xax]
        #py=pos[:,yax]  
        #pz=pos[:,zax]

        px = diskCoord[:,xax]
        py = diskCoord[:,yax] 
        pz = diskCoord[:,zax] 

        px = np.append(px, bulgeCoord[:,xax])
        py = np.append(py, bulgeCoord[:,yax])
        pz = np.append(pz, bulgeCoord[:,zax])
       
        
        #xcm = px.sum()/len(px)
        #ycm = py.sum()/len(py)
        #zcm = pz.sum()/len(pz)
        xcm = np.median(px)
        ycm = np.median(py)
        zcm = np.median(pz)

        print( xcm, ycm, zcm, len(px), len(py), len(pz))

        px = px - xcm
        py = py - ycm
        pz = pz - zcm

        #xnew = px.sum()/len(px)
        #ynew = py.sum()/len(py)
        #znew = pz.sum()/len(pz)
        
        xnew = np.median(px)
        ynew = np.median(py)
        znew = np.median(pz) 
        print( xnew, ynew, znew, len(px), len(py), len(pz))

        #pmax=0.0
        #pmin=0.0
        #if ((pmax==0.0) & (pmin==0.0)):
        #        pmin=pz[pz>-np.inf].min()
        #        pmax=pz.max()
        #else:
        #        pz[pz<pmin]=pmin
        #        pz[pz>pmax]=pmax

        #print( "min/max of log10(histogram) = ", pmin, pmax)
        #vx=vel[:,xax]
        #vy=vel[:,yax]  
        #vz=vel[:,zax]

        
        vx = diskVels[:,xax]
        vy = diskVels[:,yax]
        vz = diskVels[:,zax]

        mass = myDiskMass 
     
        r=np.sqrt(px**2. + py**2.) 
        theta=np.arctan2(py,px)

        #ur=vx*np.cos(theta)+vy*np.sin(theta)
        #vt=-vx*np.sin(theta)+vy*np.cos(theta)

      
        #kindex=(r>6.4) & (r<6.6) & (theta>2.45) & (theta<2.5)
         

        h, x, y = np.histogram2d(r,theta,bins=[BINS_r,BINS_theta],range=[[rmin,rmax],[-np.pi,np.pi]])

        #divide by area to get surface density
        for i in range(0,BINS_r):
                h[i,:]/=rmid[i]*dr*dtheta

        #fit the axisymmetric surface density            
        meanh=np.zeros(BINS_r)
        for i in range(0,BINS_r):
                meanh[i]=h[i,:].mean()
        z=np.polyfit(rmid, np.log(meanh), 1)
        Rs=-1/z[0]
        p = np.poly1d(z)
        print( "Rs = ", Rs, mass )

        #calculate residuals
        for i in range(0,BINS_r):
                #h[i,:]=(h[i,:] - np.exp(p(rmid[i]))) / np.exp(p(rmid[i]))
                h[i,:]=(h[i,:] - h[i,:].mean()) / (h[i,:].mean())
        
        #print( mass.dtype)
        #print( "xxxxxx", mass.shape, px.shape)
        #construct mass-weighted histograms
        #Z,x,y=np.histogram2d(px/Rs,py/Rs, range=[[-lengthX,lengthX],[-lengthY,lengthY]], weights=mass, bins=BINS, normed=True)
        Z,x,y=np.histogram2d(px/Rs,py/Rs, range=[[-lengthX/Rs,lengthX/Rs],[-lengthY/Rs,lengthY/Rs]], bins=BINS, normed=True)
      
        Z=np.log10(Z)
        
   
        Zmin=Z[Z>-np.inf].min()
        Zmax=Z[Z<np.inf].max()
        if ((Zmax==0.0) & (Zmin==0.0)):
                Zmin=Z[Z>-np.inf].min()
                Zmax=Z.max()
        else:
                Z[Z<Zmin]=Zmin
                Z[Z>Zmax]=Zmax

        print( "min/max of log10(histogram) = ", Zmin, Zmax)



        fig = plt.figure(1, figsize=(25.0,25.0))
        
        #left plot
        #ax = fig.add_subplot(1,2,1,title=tname+"  t="+str(round(head.time*co.UnitTime_in_Gyr*1000.0,1))+"Myr")
        
        ax = fig.add_subplot(1,2,1) #,title=tname+"  t="+str(round(myTime*co.UnitTime_in_Gyr*1000.0,1))+"Myr")
        im=ax.imshow(Z.T, vmin=Zmin, vmax=Zmax, origin='lower',interpolation='nearest', extent=[-lengthX/Rs,lengthX/Rs,-lengthY/Rs,lengthY/Rs], cmap=cm.get_cmap('jet'))
        ax.set_xlabel('x/Rs', fontsize=18, fontweight='bold')
        ax.set_ylabel('y/Rs',fontsize=18, fontweight='bold')
        plt.xticks(np.arange(-round(lengthX/Rs), round(lengthX/Rs), step=2), fontsize=15, fontweight='bold')
        plt.yticks(np.arange(-round(lengthY/Rs), round(lengthY/Rs), step=2), fontsize=15, fontweight='bold')
        plt.colorbar(im, shrink=0.35)        

#  I DEFINE THE PATCH
        ##rlow=6.0
        ##rhigh=7.0
        ##thetalow=0.0
        ##thetahigh=1.0
           
        ##add on top particle type 3
        ##kindex=(px/(2*rh)<1.03) & (px/(2*rh)>0.97) & (py/(2*rh)>0.97) & (py/(2*rh)<1.03)  
       # kindex=(r>8.0) & (r<8.2) & (theta>2.45) & (theta<5.5)  
       # print( "UFF", len(kindex))
       # urnew=ur[kindex]-vx0
        # vtnew=vt[kindex]-vy0
        ##im=ax.imshow(h.T, vmin=-0.75, vmax=0.75, origin='lower',interpolation='nearest',extent=[-300, 300, -300, 300], cmap=cm.get_cmap('jet'))
        # ax.scatter(urnew,vtnew,s=patch_s,c='w',facecolors='black',edgecolors='black')
        # if (head.npart[2] > 0) & (patch_flag):
        #       urnew=ur[kindex]-vx0
        #        vtnew=vt[kindex]-vy0
        #       ax.scatter(urnew,vtnew,s=patch_s,c='w',facecolors='black',edgecolors='black')
  


#       if (head.npart[2] > 0) & (mc_flag):
#                if (r.all < rhigh) & (r.all > rlow) & (theta.all < thetahigh) & (theta.all > thetalow):
#                    print( r,theta )
#                    poss=rs.read_block(filename, "POS ", parttype=2, arepo=arepo).astype('float64')
#                    pxx=poss[:xax]
#                    pyy=poss[:yay]
#                    ax.scatter(pxx/rh,pyy/rh,s=mc_s,c='w',facecolors='black',edgecolors='black')
#                    ax.axis([-lengthX,lengthX,-lengthY,lengthY])
                        #ax.scatter(ur,vt,s=mc_s,c='w',facecolors='black',edgecolors='black') 
                ##pos=rs.read_block(filename, "POS ", parttype=3, arepo=arepo).astype('float64')
                #poss=rs.read_block(filename, "POS ", parttype=3).astype('float64')
                #pxx=poss[:,xax]
                #pyy=poss[:,yax]
                #ax.scatter(pxx/rh,pyy/rh,s=mc_s,c='w',facecolors='black',edgecolors='black')
                #ax.axis([-lengthX,lengthX,-lengthY,lengthY])
        ##if counter==len(nums1)-1:
        #print( 'nnn', len(urnew))

# commento da Tenerife
        # ax.set_xlabel('u', fontsize=18)
        # ax.set_ylabel('v', fontsize=18)
# fine commento

        #print( 'cccc', pxx,pyy )
        ax2 = fig.add_subplot(1,2,2,title="residuals")
        #im=ax.imshow(h.T, vmin=-0.75, vmax=0.75, origin='lower',interpolation='nearest',extent=[0, lengthX, -np.pi,np.pi ], cmap=cm.get_cmap('jet'))
        im2=ax2.imshow(h.T, vmin=-0.75, vmax=0.75, origin='lower',interpolation='nearest',extent=[0, lengthX/Rs, -np.pi,np.pi ], cmap=cm.get_cmap('jet'))
        #add on top particle type 3
        #if (head.npart[4] > 0) & (mc_flag):
                ##pos=rs.read_block(filename, "POS ", parttype=3, arepo=arepo).astype('float64')
                #pos=rs.read_block(filename, "POS ", parttype=3).astype('float64')
                #px=pos[:,xax]
                #py=pos[:,yax]
                #r=np.sqrt(px**2. + py**2.)
                #theta=np.arctan2(py,px)
                #ax.scatter(r/Rs,theta,s=mc_s,c='w',facecolors='black',edgecolors='black')
        ax2.axis([0,lengthX/Rs,-np.pi,np.pi], thicks='')
        plt.xticks(np.arange(0, round(lengthX/Rs), step=1), fontsize=15, fontweight='bold')
        plt.yticks(np.arange(-np.pi, np.pi, step=1), fontsize=15, fontweight='bold')

        ax2.set_xlabel('R/Rs', fontsize=18, fontweight='bold')
        ax2.set_ylabel(r'$\theta$',fontsize=20, fontweight='bold')
        cbar2=plt.colorbar(im2, shrink=0.35)
        cbar2.set_label("$(\Sigma - \Sigma_{\mathrm{fit}}) / \Sigma_{\mathrm{fit}} $", fontsize=18, fontweight='bold')

        #plt.colorbar(im, shrink=0.5)
###     plt.savefig(output_dir+"uvmultipatch_8_150_"+str(num).zfill(3)+".png", dpi=100, bbox_inches='tight')
        plt.savefig(output_dir+"M4disk_new"+str(num).zfill(3)+".png", dpi=100, bbox_inches='tight')
        #plt.savefig("bar306.pdf") 
        fig.clf()




##############################################

#things to do to read hdf5 for gadget

"""

#this is for a single snapshot

myFile = h5py.File("name_of_my_simulation_snap....hdf5","r")

#if you want to look at header:

myHead = myFile.__getitem__("/Header")

#same with part types:

myPart1 = myFile.__getitem__("PartType1")

myCoord1 = np.array(myPart1.__getitem__("Coordinates"))

#same with velocity, acceleration.....


















"""




