# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 10:31:13 2015

@author: Student
"""

from __future__ import division
import pygame   #Imports methods from a file called 'pygame'
from vector import *   #Rename your solution to homework 2 as vector.py 
import random      #Needed to generate random numbers
import math
from itertools import combinations
from random import randint
import ode

author_name = "Jeremy Rodgers"

class particle:
	def __init__(self,pos, vel, mass = 1.0, size = 10, shape = "circle", color = (0,0,255)):
		self.pos = pos   #Position vector
		self.vel = vel   #Velocity vectorB
		self.m = mass    #Particle mass
		self.m_inv = 1/mass
		self.size = size #Radius for circle and width and height for rectangle
		self.color = color
		self.shape = shape  #String with the name of the shape
		self.f_net = vect2d(0,0)   #Property to keep track of net force on particle
		self.pos_old = vect2d(0,0)
		self.vel_old = vect2d(0,0)
		self.f_net_old = vect2d(0,0)
		self.k1 = vect2d(0,0)
		self.k2 = vect2d(0,0)
		self.k3 = vect2d(0,0)
		self.k4 = vect2d(0,0)
		self.L1 = vect2d(0,0)
		self.L2 = vect2d(0,0)
		self.L3 = vect2d(0,0)
		self.L4 = vect2d(0,0)
		self.drag_coeff = 0.000001   #Scaling number for drag force
		self.cor = .99     #Coefficient of restitution
		self.rod_connection={}  #NEW: Dictionary of all rod constraint connections
		self.geom = [] #Place holder for the Geom


	def add_force(self,force):
		"""Add force to the net force vector for each particle"""
		self.f_net += force


	def reset_force(self):
		"""Set net force to zero

		Call this at the start of each time step to reset net forces for each particle
		"""
		self.f_net = vect2d(0,0) 


	def get_force(self):
		"""Return f_net property"""
		return self.f_net

	def set_pos(self,new_pos):
		"""Set the position of the particle

		Make sure that new_pos is a vect2d object
		"""
		self.pos = new_pos
	
	def set_gCoord(self):
		
		self.geom.setPosition(self.pos.x, self.pos.y, 0)



class world:
	def __init__(self):
		self.particle_list = []
		self.all_particles = []
		self.background_color = (255,255,255)  #World background set to white
		self.dt = 0.1     #time step size - Fixed for now
		self.width = 800  
		self.height = 600
		self.force_que = []    #Keeps track of forces, particles, and things like spring constants, etc
		self.g = vect2d(0,1)  #Constant gravitational field constant
		self.numerical = 'rk4'
		self.running = True   #Determines if while look should continue
		self.selected = None   #Particle selected with the mouse
		self.mouse_force = 1    #Coefficient for force exerted by mouse
		self.damping = 0.5    #Damping coefficient
		self.vel_max = 50     #Velocity at which the drag force kicks in
		self.rod_color=(120,120,0)  #NEW: Used for color of connecting rods
		self.setup_world()   #Create the pygame window


	def setup_world(self):
		"""Create a pygame window"""
		pygame.init()
		self.screen = pygame.display.set_mode((self.width,self.height))
		pygame.display.set_caption('PHYS-360 Homework 9')
		self.screen.fill(self.background_color)
		self.clock = pygame.time.Clock()
		self.collision = collision_engine(self,self.particle_list)
		pygame.time.Clock()



	def set_numerical(self,method):
		"""Change the self.numerical variable to use a different update method

		Note that 'method' must be passed as a string
		"""
		self.numerical = method

	def add_particle(self,particle):
		"""Add particle to particle_list"""
		self.particle_list.append(particle)
		self.all_particles.append(particle)
		
		particle.geom = ode.GeomSphere(self.collision.space, particle.size)
		particle.geom = ode.GeomBox(self.collision.space,lengths = particle.width, particle.height,1)


	def add_rod(self,part1,part2):
		"""NEW: Adds a rod constraing between two particles with length equal to separation at initial time

		The length of the rod is set to the original distance between the two particles.  Particle class
		has new property called rod_connection.  This is a dictionary that uses the name of the particle as a key and the length of the rod 
		as the associated value.  You may want to look up python dictionary. 
		"""
		L = (part1.pos - part2.pos).mag()  #Length of rod
		part1.rod_connection[part2] = L
		part2.rod_connection[part1] = L


	def new_force(self,force,particles):
		"""Add a force to the list of forces calculated

		'particles' should be a list of all particles experiencing the force.  Some
		forces require other parameters and those should be included in the 'particles' list
		(e.g. see spring force defintion below)
		"""
		self.force_que.append([force,particles])


	def update(self):
		"""Update the positions and velocities of all particles

		Should also include boundary(), self.game_controls(), self.v_max_check(),
		and should use getattr() to call the numerical method used to solve
		the equations of motion.
		"""
		self.game_controls()
		getattr(self,self.numerical)()  #Solve equations of motion
		self.collision.find_collision()  #Check to see if object hits the walls
		self.display()


	def display(self):
		"""Draw all particles onto the surface and then flip to computer screen"""
		self.screen.fill(self.background_color)
		for i in self.particle_list:
			if i.shape == "rectangle":
				i.draw(self.screen)	
			else:
				pygame.draw.circle(self.screen, i.color, (int(i.pos.x), int(i.pos.y)), i.size, 1)

			pygame.display.flip()




	def move(self):
		"""Updates the position and velocity of particles under influence of a changing force"""
		for i in self.particle_list:
			i.vel += i.f_net*self.dt
			i.pos += i.vel*self.dt
			i.set_gCoord()


	def euler(self):
		"""Update the position and velocity using the Euler method"""
		for i in self.particle_list:
			i.vel += i.f_net*self.dt
			i.pos += i.vel*self.dt
			i.set_gCoord()

			if i.shape == "rectangle":
				i.apply_rotations(self.dt)




	def midpoint(self):
		"""Update the position and velocity using the Midpoint method"""
		for i in self.particle_list:
			velInit = i.vel			
			i.vel += i.m_inv*i.f_net*self.dt
			i.pos += (1/2)*(velInit+i.vel)*self.dt
			i.set_gCoord()


	def verlet(self):
		"""Use the Velocity Verlet method to update position and velocity"""
		for i in self.particle_list:
			i.pos += i.vel + (1/2)*(i.m_inv*i.f_net)*self.dt**2
			i.fNetOld = i.f_net
			i.set_gCoord()



	def projectile(self):
		"""Updates the position of a particle moving with a constant net force

		Make sure that other position update methods are not called at the same time
		as this one (e.g. either projectile() or move(), not both).
		"""
		for i in self.particle_list:
			i.pos += i.vel*self.dt + 1/2*self.g*self.dt**2
			i.vel += self.g*self.dt
			i.set_gCoord()



	def rk4(self):
		"""Update the position and velocity using rk4

		Will need to use particle.get_force() method and particle.set_pos() at different points to 
		return the net force on each particle and set the positions before calling self.net_force().
		Store k and L results (use capital L so it isn't confused with number one)"""
		for i in self.particle_list:
			i.pos_old = i.pos
			i.vel_old = i.vel
			i.f_net_old = i.f_net
			
			i.k1 = i.vel*self.dt
			i.L1 = (i.f_net/i.m)*self.dt			
			i.k2 = (i.vel + (1/2)*i.L1)*self.dt
			#update position
			i.pos = i.pos_old + (i.k1/2)
		#Re-calculate net forces on all particles
		self.net_force()
		
		for i in self.particle_list:
			i.L2 = (i.f_net/i.m)*self.dt			
			i.k3 = (i.vel + (1/2)*i.L2)*self.dt
			#update position
			i.pos = i.pos_old + (i.k2/2)
		#Re-calculate net forces on all particles
		self.net_force()

		for i in self.particle_list:			
			#update position
			i.L3 = (i.f_net/i.m)*self.dt			
			i.k4 = (i.vel + i.L3)*self.dt
			#update position
			i.pos = i.pos_old + i.k3
		self.net_force()
		
		for i in self.particle_list:	
			self.L4 = (i.f_net/i.m)*self.dt
			#reset to old
			i.pos = i.pos_old
			i.f_net = i.f_net_old
			i.vel = i.vel_old
			#Final rk4 Equation
			i.vel += (1/6)*i.L1 + (1/3)*i.L2 + (1/3)*i.L3 + (1/6)*i.L4
			i.pos += (1/6)*i.k1 + (1/3)*i.k2 + (1/3)*i.k3 + (1/6)*i.k4
			i.set_gCoord()

			if i.shape == "rectangle":
				i.apply_rotations(self.dt)



	def net_force(self):
		"""Find net force for all particles

		Set net force for each particle to zero at start of the time setup_world
		and then run through the force_que to find net force on each particle
		"""
		for i in self.particle_list:
			i.reset_force()   #Reset all net forces to zero
		for f in self.force_que:
			getattr(self,f[0])(f[1])  #Calls all force methods listed in force_que


	def const_grav(self,particles):
		"""Constant gravitational field

		Don't forget that the gravitational FORCE depends on the mass
		"""
		for i in particles:
			F = self.g*i.m
			i.add_force(F)


	def drag(self,particles):
		"""Apply drag force to particles

		Drag force only kicks in when particle speed exceeds self.vel_max.
		Drag force is opposite direction of velocity and should scale as
		the square of the speed
		"""
		for i in particles:
			if i.vel.mag() > self.vel_max:
				i.add_force(-i.drag_coeff*i.vel.mag2()*i.vel)

	def spring(self,particles):
		"""Spring force between two particles

		Assumes only two particles are passed.  If you have more than
		two particles you will need to midify the code.  The argument 'particles' should be a
		list containing the two particles, the spring constant, and the relaxed separation distance
		"""
		L = particles[0].pos - particles[1].pos
		F = -particles[2]*(L.mag() - particles[3])*L.norm()		
		#F = -particles[2]*(L.mag() - particles[3]) + self.damping*(particles[0].vel.x - particles[1].vel.x)*L.norm()
		particles[0].add_force(F)
		particles[1].add_force(-F)


	def mouse_pull(self,particle):
		"""Force applied by selecting particle with mouse"""
		(pick_x,pick_y) = pygame.mouse.get_pos()
		mouse_pos = vect2d(pick_x,pick_y)
		dx = self.selected.pos - mouse_pos
		F = -self.mouse_force*dx - self.selected.vel*self.damping
		self.selected.add_force(F)



	def game_controls(self):
		"""Handles all mouse and keyboard events

		Call this in the update() method
		"""
		for event in pygame.event.get(): 
			if event.type == pygame.QUIT: #If red 'x' is clicked
				self.running = False
			if event.type == pygame.MOUSEBUTTONDOWN:
				(pick_x,pick_y) = pygame.mouse.get_pos()   #Get mouse position on screen
				picked = vect2d(pick_x,pick_y)  #Turn mouse position into a vector
				for i in self.particle_list:
					dist = i.pos - picked   #How far the mouse is from the center of a particle
					if dist.mag() < i.size:  #If mouse click is inside circle
						self.selected = i
						self.new_force('mouse_pull',[self.selected])  #Add mouse force to force_que
						self.selected.color = (255,0,0)   #Change color
			if event.type == pygame.MOUSEBUTTONUP and self.selected != None:
				self.selected.color = (0,0,255)  #change color back to blue
				self.force_que.remove(['mouse_pull',[self.selected]])   #Remove force from force_que
				self.selected = None    #No particles are selected


class collision_engine():
	def __init__(self,world, particle_list=[]):
		self.world = world  #World instance used by this collision engine
		self.particle_list = particle_list   #List of all particles in the world
		self.has_collided = []    #List of particles that have collided with each other
		self.cor = 0.99   #Coefficient of restitution for hitting walls
		self.space = ode.Space()


	def find_collision(self):
		"""EDIT: Runs through all particles in particle_list to look for overlapping particles

		Needs to check and see if particle has a rod connection and if it does, call check_rod_constraint
		"""
		self.has_collided = []   #Reset collisions
		self.boundary()    #See if any particles collide with the wall
		for p in combinations(self.particle_list,2): 
			contact = ode.collide(p[0].geom, p[1].geom) 
			if contact != []: 
				r_c = contact[0].getContactGeomParams()[0] 
				n_hat = contact[0].getContactGeomParams()[1]
				dx = contact[0].getContactGeomParams()[2] 
				self.resolve_collision_ode(p[0],p[1],r_c,n_hat,dx)


	def bounding_sphere(self,particle1,particle2):
		"""Check for overlap using bounding spheres"""
		if particle1.shape == "rectangle":
			rRel = particle1.pos - particle2.pos
			
			xBodyCoord = rRel*particle1.xBodyCoord
			yBodyCoord = rRel*particle1.yBodyCoord

			if (xBodyCoord <= (particle1.width/2)) and (xBodyCoord >= (-particle1.width/2)) and (yBodyCoord <= (particle1.height/2)) and (yBodyCoord >= (-particle1.height/2)):
				particle1.collisionCoord.x = xBodyCoord
				particle1.collisionCoord.y = yBodyCoord
				self.has_collided.append([particle1, particle2])
		
		elif (particle1.shape == "circle") and (particle2.shape == "circle"):
			if abs((particle1.pos - particle2.pos).mag()) < particle1.size + particle2.size:
				self.has_collided.append([particle1, particle2])



	def resolve_collision(self):
		"""Changes velocities of point particles that have collided"""
		for i in self.has_collided:
			
			if i[0].shape == "rectangle":
				
				velPart1 = ((1/i[0].m)*(1 + self.cor)*(i[1].vel - i[0].vel))/((1/i[1].m) + (1/i[0].m))
				velPart2 = -velPart1
				
				i[0].vel += velPart1
				i[1].vel += velPart2

				
				i[0].collisionCoord = i[0].collisionCoord.x*i[0].pos.norm() + i[0].collisionCoord.y*i[0].pos.norm()
				hat = i[0].collisionCoord.x*i[0].xBodyCoord.norm() + i[0].collisionCoord.y*i[0].yBodyCoord.norm()

				Vclose = i[1].vel - i[0].vel

				collisionCoordPerp = vect2d(0,0)
				collisionCoordPerp.x = i[0].collisionCoord.x*(math.cos(90*(math.pi/2))) - i[0].collisionCoord.x*(math.sin(90*(math.pi/2)))
				collisionCoordPerp.y = i[0].collisionCoord.y*(math.sin(90*(math.pi/2))) + i[0].collisionCoord.y*(math.cos(90*(math.pi/2)))

				I = (1/2)*i[0].m*((i[0].omega**2) + collisionCoordPerp.mag2())
				print '{}'.format(I)
				
				J = ((1 + self.cor)*Vclose.mag())/(((1/i[1].m) + (1/i[0].m)) + ((collisionCoordPerp*hat)/I))
				print'{}'.format(J)
				print '{}'.format((-collisionCoordPerp*hat*J)/I)

				i[0].omega = (-collisionCoordPerp*hat*J)/I

			elif (i[0].shape == "circle") and (i[1].shape == "circle"):
				velPart1 = ((1/i[0].m)*(1 + self.cor)*(i[1].vel - i[0].vel))/((1/i[1].m) + (1/i[0].m))
				velPart2 = -velPart1
				
				i[0].vel += velPart1
				i[1].vel += velPart2

			self.resolve_interpenetration(i[0], i[1], self.collision_normal(i[0], i[1]))

	def resolve_collision_ode(self):

		print'Finish'


	def collision_normal(self,particle1,particle2):
		"""Returns the collision normal vector"""
		return -(particle1.pos-particle2.pos).norm()



	def resolve_interpenetration(self,particle1,particle2,hat):
		"""Move objects that overlap far enough apart that they don't continue to collide"""
		dx = particle1.size + particle2.size - (particle1.pos - particle2.pos).mag()
		
		particle1.pos += ((particle1.m/(particle1.m + particle2.m))*dx)*hat
		particle2.pos += ((particle2.m/(particle1.m + particle2.m))*dx)*hat


	def boundary(self):
		"""Checks to see if part of a particle leaves the screen.  Should shift particle back on to screen and reverse component of velocity

		This method has been moved from the world class.  It makes more sense to have it be part of the collision engine
		"""

		for i in self.particle_list:
			
			if i.shape == "rectangle":
				if i.pos.x + i.width/2 >= self.world.width:			
					i.vel.x *= -self.cor
					i.pos.x += i.vel.x
				if i.pos.y + i.height/2 >= self.world.height:
					i.vel.y *= -self.cor
					i.pos.y += i.vel.y
				if i.pos.x - i.width/2 <= 0:
					i.vel.x *= -self.cor
					i.pos.x += i.vel.x
				if i.pos.y - i.height/2 <= 0:
					i.vel.y *= -self.cor
					i.pos.y += i.vel.y

			elif i.shape == "circle":
				if i.pos.x + i.size >= self.world.width:			
					i.vel.x *= -self.cor
					i.pos.x += i.vel.x
				if i.pos.y + i.size >= self.world.height:
					i.vel.y *= -self.cor
					i.pos.y += i.vel.y
				if i.pos.x - i.size <= 0:
					i.vel.x *= -self.cor
					i.pos.x += i.vel.x
				if i.pos.y - i.size <= 0:
					i.vel.y *= -self.cor
					i.pos.y += i.vel.y


	def check_rod_constraint(self,part0, part1,L):
		"""NEW: Check to see if particles are too close or too far apart.  Call fix_rod_constraint if Needed

		L is the length of the rod.  To get the length L you will need to use part0.rod_connection[part1].  This uses the rod_connection dictionary to look
		up the associated value, which is the rod length L.
		You will need to calculate the vector along the rod (hat) and
		the current distance between the particles (dist) to pass to fix_rod_constraint.
		"""
		if L != ((part0.pos - part1.pos).mag() + L):
			rhat = (part0.pos - part1.pos).norm()
			dist = (part0.pos - part1.pos).mag()
			self.fix_rod_constraint(part0, part1, L, rhat, dist)
		

	def fix_rod_constraint(self,part0,part1,L,hat,dist):
		"""NEW: Move particle and change velocity so rod constraints are maintained

		L is length of rod, hat is the vector pointing along the rod, and dist is the
		current distance between the particles.  You can use the same impulse equations
		used during collisions to change the velocities.  To move the particles use the 
		interpentration resolution formula.
		"""
		c = part0.cor * part1.cor
		dvel = (part1.vel - part0.vel)*hat
		
		part0.vel += ((part0.m_inv*(1+c)*dvel)/(part1.m_inv + part0.m_inv))*hat
		part1.vel -= ((part1.m_inv*(1+c)*dvel)/(part1.m_inv + part0.m_inv))*hat
		
		dx = dist - L
		
		part0.pos += -((part0.m/(part0.m_inv+part1.m_inv))*dx)*hat
		part1.pos -= -((part1.m/(part0.m_inv+part1.m_inv))*dx)*hat
		
		pygame.draw.aaline(self.world.screen, self.world.rod_color, (part0.pos.x, part0.pos.y), (part1.pos.x, part1.pos.y), True)
		pygame.display.flip()
				

#Work Here	
class rectangle(particle):
	def __init__(self, pos, vel, mass, width, height, theta = 0, omega = 0, torque = 0, shape = "rectangle", color = (0, 255, 120)):
		"""net_torque"""
		"""net_torque_old"""
		size = vect2d(width/2, height/2).mag()

		particle.__init__(self, pos, vel, mass, size, shape, color)
		
		self.width = width
		self.height = height
		self.theta = theta
		self.omega = omega
		self.torque = torque
		self.xBodyCoord = vect2d(1,0)
		self.yBodyCoord = vect2d(0,1)
		self.collisionCoord = vect2d(0,0)

		
	def draw(self, screen):
		corners = self.find_corners()
		
		pygame.draw.aalines(screen, self.color, True, [(corners[0].x, corners[0].y), (corners[1].x, corners[1].y), (corners[2].x, corners[2].y), (corners[3].x, corners[3].y)])

	def net_mass(self, mass):
		self.mass = mass
		self.mass.m_inv = 1/mass
		
	def find_corners(self):

		corner1 = vect2d((self.pos.x + ((-width/2)*self.xBodyCoord.x + (-height/2)*self.xBodyCoord.y)), (self.pos.y + ((-width/2)*self.yBodyCoord.x + (-height/2)*self.yBodyCoord.y)))
		corner2 = vect2d((self.pos.x + ((-width/2)*self.xBodyCoord.x + (height/2)*self.xBodyCoord.y)), (self.pos.y + ((-width/2)*self.yBodyCoord.x + (height/2)*self.yBodyCoord.y)))
		corner3 = vect2d((self.pos.x + ((width/2)*self.xBodyCoord.x + (height/2)*self.xBodyCoord.y)), (self.pos.y + ((width/2)*self.yBodyCoord.x + (height/2)*self.yBodyCoord.y)))
		corner4 = vect2d((self.pos.x + ((width/2)*self.xBodyCoord.x + (-height/2)*self.xBodyCoord.y)), (self.pos.y + ((width/2)*self.yBodyCoord.x + (-height/2)*self.yBodyCoord.y)))
		
		cornerArray = [corner1, corner2, corner3, corner4]

		return cornerArray
		
	def old_torque(self, torque):
		torque = self.torque
		return torque
		
	def rest_torque(self):

		return self.torque
		
	def get_torque(self):
		return self.torque

	def apply_rotations(self, dt):

		self.theta += self.omega*dt

		self.xBodyCoord = vect2d(math.cos(self.theta*(math.pi/2)), -math.sin(self.theta*(math.pi/2)))
		self.yBodyCoord = vect2d(math.sin(self.theta*(math.pi/2)), math.cos(self.theta*(math.pi/2)))

	def set_gCoord(self):
		
		self.geom.setPosition((self.pos.x,self.pos.y,0)) 
		self.geom.setRotation([cos(self.theta), -sin(self.theta), 0, sin(self.theta), cos(self.theta),0,0,0,1])





	
	"""Update existing numerical fuctions to check object shape and preform
	additional calculations"""

if __name__ == '__main__': 
	#The following is just to test your code out and make sure it works
	earth = world()
	part_list = []

	earth.set_numerical('rk4')  #Change to the Velocity Verlet method

	pos = vect2d(100,100)
	vel = vect2d(3,3)
	width = 100
	height = 100
	mass = 1.0
	rect = rectangle(pos,vel,mass,width,height) #create new rectangle
	part_list.append(rect)
	earth.add_particle(rect)


	num_part = 2  #Set number of particles to add to world
	part_list = []
	for i in range(num_part):
		r = vect2d(randint(0,700),randint(0,600))
		v = vect2d(randint(-5,5),randint(-5,5))
		part = particle(r,v)   #Create new particle
		part_list.append(part)    #Add to list of particles that have been randomly placed
		earth.add_particle(part)  #Add particle to the world
	earth.new_force('const_grav',part_list)   #Particles experience gravity
	earth.new_force('drag',part_list)   #particles also experience drag



	###
	#While loop
	###

	while earth.running:   
		time_passed = earth.clock.tick(120)   #Add pygame.time.Clock() to world setup.  This will allow you to set fps
		earth.update()
	#Quit The Window
	pygame.quit()
