#Aidan Walk
#Astr 260-001
#HW09 Problem 2 part a
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

def analytic_projectile_nodrag(t,
                               initial_position=None,
                               initial_velocity=None):
    #trick to apply to x and y component simultaneously
    times = t[:,np.newaxis]
    return initial_position +\
           initial_velocity*times + \
           g*0.5*times**2
           
def accelerationOfBall(mass, t=None):
    '''equation for the acceleration of the ball'''
    acceleration = forceOfDrag()/m + g
    return acceleration
    
def forceOfDrag():
    '''
    '''
    return 0
         
def rk4(xnaught = None,
        vnaught = None,
        theta = None,
        times = None,
        h = None,
        deriv = None,
        deriv_params = None):
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
        k1v = h*deriv(v1, t=t)
        #compute k2
        v2 = currentV + 1/2*k1v
        k2v = h*deriv(v2, t=t+h/2)            
        #compute k3
        v3 = currentV + 1/2*k2v
        k3v = h*deriv(v3, t=t+h/2)
        #compute k4
        v4 = currentV+k3v
        k4v = h*deriv(v4, t=t+h)
        
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
    
    angle = 45*np.pi/180 #degrees
    
    speed = 10 #m/s
    x0 = np.array([0,2])
    v0 = speed*np.array([np.cos(angle), np.sin(angle)])
    
    
    #compute analytic solution with no drag
    analyticNoDrag = analytic_projectile_nodrag(times,
                                           initial_position=x0,
                                           initial_velocity=v0)
    #compute position via RK4                                   
    RK4Pos, RK4Vel = rk4(xnaught=x0,
              vnaught=v0,
              times = times,
              h = timeStep,
              theta = angle,
              deriv = accelerationOfBall,
              deriv_params = {'Vnaught':v0, 'a':g, 't':times})
    
    
    #plot analytic solution 
    plt.plot(analyticNoDrag[:,0], analyticNoDrag[:,1], #where positions[0] is x, [1] is y
             label='Analytic',
             color='#ff0000', linewidth=2) 
    #Plot RK4 solution         
    plt.plot(RK4Pos[:,0], RK4Pos[:,1],
             label='RK4',
             color = '#000000', linewidth=2, linestyle=':')
    
    #make plot look pretty
    plt.title('Projectile Motion - No Drag')
    plt.xlabel('x-position')
    plt.xlim(0, 15)
    plt.ylabel('y-position')
    plt.ylim(0, 5)
    legend = plt.legend(loc='upper right')
    legend.get_frame().set_facecolor('#f8f9f9')
    #save plot
    plotFName='AidanWalk_HW09_2a_Plot'
    plt.savefig(plotFName +'.png')
    print('Plot saved to:', plotFName)
    
    
    
    
    
    
    
    
    
