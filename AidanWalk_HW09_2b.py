#Aidan Walk
#Astr 260-001
#HW09 Problem 2 part b
#31 March 2021, 17:00

import numpy as np
import matplotlib.pyplot as plt
import ipdb

#Global Variables
m = 0.145 #kg, mass of baseball
g = np.array([0, -9.8]) #m/s^2, accel of grav at surface
rho = 1.225 #kg/m^3, density of air at STP
Cd = 0.47 #[unitless], drag coefficient of sphere
A = 0.00414 #m^2, cross-sectional area of baseball
mphToms = 0.44704 #meters/second = 1 mph
           
def accelerationOfBall(velocity, drag):
    '''equation for the acceleration of the ball
       velocity = np.array([x-component, y-component]) (only used if drag==1)
       drag = bool; true for air resistance to apply, else no drag
       returns acceleration of the ball'''
    if drag==True:
        acceleration = forceOfDrag(velocity)/m + g
    else:
        acceleration = g
    return acceleration
    
def forceOfDrag(v):
    '''Describes the forces on the ball due to air resistance
       take in v, array of x and y components of velocity, 
       returns Force due to drag as an array'''
    Vmag = v / (np.sqrt(v[0]**2+v[1]**2))
    Fdrag = -1/2*rho* np.abs(v)**2 *Cd*A*Vmag
    return Fdrag
         
def rk4(xnaught = None,
        vnaught = None,
        times = None,
        h = None,
        deriv = None,
        drag = None):
    '''Runge-Kutta 4 method
       returns times, vals
       vnaught = array[x,y], x and y components of initial velocity
       '''
    #initialize arrays where data will be stored as well as vnaught, xnaught
    velcities = currentV = vnaught
    positions = currentPos = xnaught
    #compute velocity, position for each time in times
    for t in times:
        #update arrays with current information
        positions = np.vstack((positions, currentPos))
        velocities = np.vstack((velcities, currentV))
        
        #calculate velocity:
        #compute k1 
        v1 = currentV
        k1v = h*deriv(v1, drag)
        #compute k2
        v2 = currentV + 1/2*k1v
        k2v = h*deriv(v2, drag)            
        #compute k3
        v3 = currentV + 1/2*k2v
        k3v = h*deriv(v3, drag)
        #compute k4
        v4 = currentV+k3v
        k4v = h*deriv(v4, drag)
        
        #calculate position
        k1x = (currentV) *h
        k2x = (currentV + k1v*h/2) *h
        k3x = (currentV + k2v*h/2) *h
        k4x = (currentV + k3v*h) *h       
        
        #update current velocity and position
        currentV = currentV + 1/6 * (k1v + 2*k2v + 2*k3v + k4v)
        currentPos = currentPos + 1/6 * (k1x + 2*k2x + 2*k3x + k4x)
        
    return positions, velcities

if __name__ == "__main__":
    maxTime = 10
    timeStep = 0.01
    times = np.arange(0, maxTime, timeStep)
    speed = 85 *mphToms #m/s
    x0 = np.array([0,2])
    
    #compute RK4 for each angle, with and without drag
    angles = [60, 45, 30]
    drags = [True, False]
    for angle in angles:
        angleRad = angle*np.pi/180 #convert degrees to radians
        v0 = speed*np.array([np.cos(angleRad), np.sin(angleRad)]) 
        
        for drag in drags:
            #compute position via RK4
            RK4Pos, RK4Vel = rk4(xnaught=x0,
                                 vnaught=v0,
                                 times = times,
                                 h = timeStep,
                                 deriv = accelerationOfBall,
                                 drag = drag)
                                 
            #assign plot characteristics based on angle and drag
            if angle == 30: color = '#347502'
            if angle == 45: color = '#026075'
            if angle == 60: color = '#a10000'
            width=2
            if drag==False: #If there is no drag, make lines dotted and black, exclude labels
                style = ':' 
                color = 'k'  
                labels = '' 
                width = 1
            else: #if there is drag, assign label and linestyle
                style = '-'
                labels = str(angle)+'\u00b0'
            #create single label for the no drag lines
            if (drag==False and angle==30): labels = 'No Drag' 
            
            #plot the data
            plt.plot(RK4Pos[:,0], RK4Pos[:,1],
                     label=labels,
                     color = color, 
                     linewidth=width, 
                     linestyle=style)
    
    
    #make plot look pretty
    plt.title('Projectile Motion With Drag')
    plt.xlabel('x-position')
    plt.xlim(0, 150)
    plt.ylabel('y-position')
    plt.ylim(0, 60)
    legend = plt.legend(loc='upper right')
    legend.get_frame().set_facecolor('#f8f9f9')
    
    #save plot
    plotFName='AidanWalk_HW09_2b_Plot'
    plt.savefig(plotFName +'.png')
    print('Plot saved to:', plotFName)