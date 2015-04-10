
"""
Created on Fri Jan 23  2015

@author: Todd Zimmerman

Created by Todd Zimmerman for PHYS-360 Spring 2015 at UW-Stout

Used for Stem expo project
"""

import pygame   #Imports methods from a file called 'pygame'
from vector import *   #Rename your solution to homework 2 as vector.py and make sure it is the same directory as this file
import random      #Needed to generate random numbers 

author_name = "Aaron Aumann"


class particle:
    def __init__(self,pos, vel, mass = 1, size = 10, shape = "circle", color = (0,0,255)):
        """You can set default values for variables.  If you don't specify mass, size, etc, then the default
        values are used.
        """
        self.pos = pos   #Position vector
        self.vel = vel   #Velocity vectorB
        self.m = mass    #Particle mass
        self.size = size #Radius for circle and width and height for rectangle
        self.color = color
        self.shape = shape  #String with the name of the shape


class world:
    def __init__(self):
        """World will do all of the heavy lifting for us, including
        keeping track of all objects, updating position, and updating the display.
        """
        self.particle_list = []  #a list of all particle in the world
        self.background_color = (255,255,255)  #World background set to white
        self.dt = 0.1   #Set the time step
        self.setup_world()   #Create the pygame window

    def setup_world(self):
        """Create a pygame window

        This method should get pygame up and running and create a window"""
        #Add your code here
        pygame.init()#Initilize pygame
        self.width = 500
        self.height = 500
        (width, height) = (500,500)
        self.screen = pygame.display.set_mode((width,height))#create a surface in memory
        #self.screenXZ = pygame.display.set_mode((width, height))
        #self.screenYZ = pygame.display.set_mode((width, height))
        pygame.display.set_caption('Program homework 3')
        self.screen.fill(self.background_color)
        

    def add_particle(self,particle):
        """Add particle to particle_list

        Use the append function to add a particle to particle_list"""
        self.particle_list.append(particle)


    def update(self):
        """Update the positions (and eventually velocities) of all particles

        Should update position of all particles in the world.  This method will ultimately
        be used to find net force an update the velocities, but for now we will
        assume constant velocities"""
        #Add your code here
        for particle in self.particle_list:
            particle.pos +=  particle.vel
            
            #Keep the particles in the bounds of the screen, simply reverses the
            #direction of the particles once the hit edge of screen.  Will
            #be replaced later with actual boundry function
            if particle.pos.x + particle.size >= self.width:
                particle.pos.x = self.width - particle.size
                particle.vel.x *= -1
            if particle.pos.x - particle.size <= 0:
                particle.pos.x = 0 + particle.size
                particle.vel.x *= -1
            if particle.pos.y + particle.size >= self.height:
                particle.pos.y = self.height - particle.size
                particle.vel.y *= -1
            if particle.pos.y - particle.size <= 0:
                particle.pos.y = 0 + particle.size
                particle.vel.y *= -1
            if particle.pos.z + particle.size >= self.height:
                particle.pos.z = self.height - particle.size
                particle.vel.z *= -1
            if particle.pos.z - particle.size <= 0:
                particle.pos.z = 0 + particle.size
                particle.vel.z *= -1
        

    def display(self):
        """Draw all particles onto the surface and then flip to computer screen
        
        Clear the screen and draw the particles in their new locations on the surface,
        then update the computer screen
        """
        #Add your code here
        self.screen.fill(self.background_color)#Clear the screen
        #self.screenXZ.fill(self.background_color)
        
        for particle in self.particle_list:#Iterate through each particle and draw it
            pygame.draw.circle(self.screen, particle.color, (int(particle.pos.x/2), int(particle.pos.y/2)), particle.size, 1)
            pygame.draw.circle(self.screen, particle.color, (int(particle.pos.x/2) +250, int(particle.pos.z/2)), particle.size, 1)
            pygame.draw.circle(self.screen, particle.color, (int(particle.pos.y/2), int(particle.pos.z/2) + 250), particle.size, 1)
            pygame.draw.aalines(self.screen, (255,0,0), True, [(0,0), (250,0), (250,250), (0,250)], True)
            pygame.draw.aalines(self.screen, (255,0,0), True, [(0,250), (250,250), (250,500), (0,500)], True)
            pygame.draw.aalines(self.screen, (255,0,0), True, [(250,0), (500,0), (500,250), (250,250)], True)
        

        pygame.display.flip()#Display the screen


if __name__ == '__main__': 
    #You will need to create an instance of the world and populate it with some particles
    World = world()
    #Add 25 particles to the world
    numberOfParticles = 1
    for n in range(numberOfParticles):
        #Generate random position and velocity vectors
        position = vect3d(random.randint(0,500),random.randint(0, 500), random.randint(0,500))
        velocity = vect3d(0.05*random.randint(-1,1), 0.05*random.randint(-1, 1), 0.05*random.randint(-1,-1))
        #Create and add particle to world using random vectors
        lParticle = particle(position, velocity);
        World.add_particle(lParticle)

    running = True  

    ###
    #While loop
    ###

    while running:   

        for event in pygame.event.get():   #This will eventually be moved into the world class but leave it alone for now
            if event.type == pygame.QUIT: 
                running = False

        #Insert your code here to update the world and then display the particles on the screen
        World.update()
        World.display()
        
    #Close out pygame
    pygame.quit()