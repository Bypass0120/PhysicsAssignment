
"""
Created on Mon, March 9 2015
@author: Todd Zimmerman
Created by Todd Zimmerman for PHYS-360 Spring 2015 at UW-Stout
Template for introducing numerical methods.  Implements euler, midpoint, and the velocity verlet methods
"""

from __future__ import division
import pygame   #Imports methods from a file called 'pygame'
from vector import *   #Rename your solution to homework 2 as vector.py 
import random      #Needed to generate random numbers
from itertools import combinations
from random import randint
import math
import ode

author_name = "Aaron Aumann"

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
		self.cor = .97     #Coefficient of restitution
		self.rod_connection={}  #NEW: Dictionary of all rod constraint connections
		self.geom = []


	def draw(self,screen, parent_pos=vect2d(0,0)):
		pygame.draw.circle(screen, self.color, (int(self.pos.x+parent_pos.x), int(self.pos.y+parent_pos.y)), self.size, 1)

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

	def set_gcoord(self):
		self.geom.setPosition((self.pos.x, self.pos.y, 0))

class rectangle(particle):
	def __init__(self, pos, vel, width=20, height=20, theta=0, omega=0, mass=1.0, shape="rectangle", color=(255, 0, 255)):
		size = vect2d(width/2, height/2).mag()
		particle.__init__(self, pos, vel, mass, size, shape, color)
		self.width = width
		self.height = height
		self.theta = theta
		self.omega = omega
		self.local_x_axis = vect2d(math.cos(theta), math.sin(theta))
		self.local_y_axis = vect2d(-math.sin(theta), math.cos(theta))
		self.torque = vect2d(0,0)
		self.torque_old = vect2d(0,0)
		self.test = False;

	def draw(self, screen):
		pygame.draw.aalines(screen, self.color, True, self.find_corners(), 1)
		#pygame.draw.circle(screen, self.color, (int(self.pos.x), int(self.pos.y)), (int(self.size)), 1)
		#Draw the local x and y axis
		#pygame.draw.aaline(screen, (255,0,0), [self.pos.x, self.pos.y], [self.pos.x + (self.width/2 * self.local_x_axis).x, self.pos.y + (self.width/2 * self.local_x_axis).y], True)
		#pygame.draw.aaline(screen, (0,0,255), [self.pos.x, self.pos.y], [self.pos.x + (self.height/2 * self.local_y_axis).x, self.pos.y + (self.height/2 * self.local_y_axis).y], True)

	def set_mass(self, mass):
		self.m = mass
		self.m_inv = 1/mass

	def find_corners(self):
		lX = self.width/2*self.local_x_axis
		lY = self.height/2*self.local_y_axis

		vCorners = [-lX -lY, lX-lY, lX+lY, -lX+lY]
		corners = [[vCorners[0].x + self.pos.x, vCorners[0].y + self.pos.y],
					[vCorners[1].x + self.pos.x, vCorners[1].y + self.pos.y],
					[vCorners[2].x + self.pos.x, vCorners[2].y + self.pos.y],
					[vCorners[3].x + self.pos.x, vCorners[3].y + self.pos.y]
		]
					
		return corners


	def add_torque(self, torque):
		self.torque += torque

	def reset_torque(self):
		print "reset"

	def get_torque(self):
		return self.torque

	def update_rotation(self, dt):
		self.theta += self.omega * dt
		#self.local_x_axis = vect2d(math.cos(self.theta * math.pi/180), -math.sin(self.theta * math.pi/180))
		#self.local_y_axis = vect2d(math.sin(self.theta * math.pi/180), math.cos(self.theta * math.pi/180))
		self.local_x_axis = vect2d(math.cos(self.theta), math.sin(self.theta))
		self.local_y_axis = vect2d(-math.sin(self.theta), math.cos(self.theta))

	def set_gcoord(self):
		self.geom.setPosition((self.pos.x, self.pos.y, 0))
		self.geom.setRotation([math.cos(self.theta), -math.sin(self.theta), 0, 
							   math.sin(self.theta), math.cos(self.theta),0,0,0,1])

class magnet(rectangle):
	def __init__(self, pos):
		rectangle.__init__(self, pos, vect2d(0,0), 100, 20)
		self.magnet_strength = 0

class coil(rectangle):
	def __init__(self, pos):
		rectangle.__init__(self, pos, vect2d(0,0), 150, 100)
		self.num_loops = 4.55; #Had to add on 0.05 to get the last point entering the coil
		self.wire = []
		self.electrons = []
		self.electron_path = []
		
		#Add visual representation for the wire not in the coil
		r = vect2d(self.pos.x - self.width/2 + 5, self.pos.y + self.height/2)
		width = 10
		height = self.height
		#Wire on the left
		self.wire.append(rectangle(r, vect2d(0,0), width, height))
		r = vect2d(self.pos.x + self.width/2 -5, self.pos.y + self.height/2)
		#Wire on the right
		self.wire.append(rectangle(r, vect2d(0,0), width, height))
		r = vect2d(self.pos.x, self.pos.y + height-5)
		#Wire on the bottom
		self.wire.append(rectangle(r, vect2d(0,0), self.width, 10))
		
		#Set up the path
		#path points for coil
		for i in range(0,(int)(self.num_loops*20)):
			self.electron_path.append(vect2d(i*(self.width-5)/(20*self.num_loops),
					    ((self.height-10)/2) * -math.sin(i*2*math.pi/20)))
		#add path points in the bottom wires
		for i in range(0, 5):#Add 5 points along the wire on the right
			self.electron_path.append(vect2d(self.width-10, 15 + i*(self.wire[1].height-15)/5))
		for i in range(0,7):#Add 7 points along the wire on the bottom
			self.electron_path.append(vect2d(self.wire[2].width-10 - i*(self.wire[2].width+10)/7,
									 		 self.wire[0].height-5))
		for i in range(0,5):#Add 5 points along the wire on the left
			self.electron_path.append(vect2d(0, 15+i*(self.wire[0].height-15)/5))

		#Add electrons to the system
		num_electrons = 100
		for i in range(0, num_electrons):
			r = self.electron_path[(int)(i*len(self.electron_path)/num_electrons)]
			r.x -= self.width/2 - 5
			part = particle(r, vect2d(0,0), 1, 5)
			self.electrons.append(part)


	def draw(self, screen):
		self.update()
		rectangle.draw(self, screen)
		for i in self.electrons:
			i.draw(screen, self.pos)
		for i in self.wire:
			i.draw(screen)

	def update(self):
		r = vect2d(self.pos.x - self.width/2 + 5, self.pos.y + self.height/2)
		width = 10
		height = self.height
		self.wire[0].pos = r
		r = vect2d(self.pos.x + self.width/2 -5, self.pos.y + self.height/2)
		self.wire[1].pos = r
		r = vect2d(self.pos.x, self.pos.y + height-5)
		self.wire[2].pos = r



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
		pygame.display.set_caption('PHYS-360 Homework 7')
		self.screen.fill(self.background_color)
		self.clock = pygame.time.Clock()
		self.collision = collision_engine(self,self.particle_list)   #This creates an instance of the collision engine



	def set_numerical(self,method):
		"""Change the self.numerical variable to use a different update method
		Note that 'method' must be passed as a string
		"""
		self.numerical = method

	def add_particle(self,particle):
		"""Add particle to particle_list"""
		if particle.shape == "circle":
			particle.geom = ode.GeomSphere(self.collision.space, particle.size)
		elif particle.shape == "rectangle":
			particle.geom = ode.GeomBox(self.collision.space, lengths=(particle.width, particle.height, 100))
			
		self.particle_list.append(particle)
		self.all_particles.append(particle)


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
		self.net_force()
		self.game_controls()
		getattr(self,self.numerical)()  #Solve equations of motion
		self.collision.find_collision()  #Check to see if object hits the walls
		self.display()


	def display(self):
		"""Draw all particles onto the surface and then flip to computer screen"""
		self.screen.fill(self.background_color)#Clear the screen
		
		for particle in self.particle_list:#Iterate through each particle and draw it
			particle.draw(self.screen)

		pygame.display.flip()#Display the screen

	def move(self):
		"""Updates the position and velocity of particles under influence of a changing force"""


	def euler(self):
		"""Update the position and velocity using the Euler method"""
		for i in self.particle_list:
			i.vel += i.m_inv * i.f_net * self.dt
			i.pos += i.vel * self.dt
			i.set_gcoord()
			#If this is a rectangle, update angular velocities as well
			if i.shape == "rectangle":
				i.update_rotation(self.dt)

	def midpoint(self):
		"""Update the position and velocity using the Midpoint method"""
		for i in self.particle_list:
			i.pos += 0.5 * (i.vel + (i.vel + (i.f_net/i.m) * self.dt)) * self.dt
			i.vel += i.f_net/i.m * self.dt
			i.set_gcoord()


	def verlet(self):
		"""Use the Velocity Verlet method to update position and velocity"""
		for i in self.particle_list:
			#Find next position
			i.pos += i.vel * self.dt + 0.5 * i.m_inv * i.f_net * self.dt**2
			#Use next position to find next force
			i.f_net_old = i.f_net
			i.set_gcoord()
			
		self.net_force()
		
		for i in self.particle_list:
			#Find next velocity using old force and next force
			i.vel += self.dt * 0.5 * i.m_inv *(i.f_net + i.f_net_old)
			#Save next force as old force


	def projectile(self):
		"""Updates the position of a particle moving with a constant net force
		Make sure that other position update methods are not called at the same time
		as this one (e.g. either projectile() or move(), not both).
		"""
		for i in self.particle_list:
			i.pos += i.vel*self.dt + 1/2*self.g*self.dt**2
			i.vel += self.g*self.dt
			i.set_gcoord()


	def rk4(self):
		"""Update the position and velocity using rk4
		Will need to use particle.get_force() method and particle.set_pos() at different points to 
		return the net force on each particle and set the positions before calling self.net_force().
		Store k and L results (use capital L so it isn't confused with number one)"""
		for i in self.particle_list:
			#Save old variables first
			i.pos_old = i.pos
			i.vel_old = i.vel
			i.f_net_old = i.f_net
			#Find k1
			i.k1 = i.vel * self.dt
			#Find L1
			i.L1 = i.m_inv * i.f_net * self.dt
			#Find k2
			i.k2 = (i.vel + 0.5 * i.L1) * self.dt
			#Move particle forwarad
			i.pos += i.k1/2
			
		#Recalc all net forces
		self.net_force()
		
		for i in self.particle_list:
			#Find L2
			i.L2 = i.m_inv * i.f_net * self.dt			
			#Find k3
			i.k3 = (i.vel + i.L2/2) * self.dt
			#Move particle forward
			i.pos = i.pos_old + i.k2/2
		
		#Recalc all net forces
		self.net_force()
		
		for i in self.particle_list:
			#Find L3
			i.L3 = i.m_inv * i.f_net * self.dt
			#Find k4
			i.k4 = (i.vel + i.L3)* self.dt
			#Move particle forward
			i.pos = i.pos_old + i.k3
		
		#Recalc all forces one more time
		self.net_force()
		
		for i in self.particle_list:
			#Find L4
			i.L4 = i.m_inv * i.f_net * self.dt
			#Reset particle to original position
			i.pos = i.pos_old
			i.vel = i.vel_old
			i.f_net = i.f_net_old
			#Find new position and velocity
			i.pos += i.k1/6 + i.k2/3 + i.k3/3 + i.k4/6
			i.vel += i.L1/6 + i.L2/3 + i.L3/3 + i.L4/6
			i.set_gcoord()


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
			i.add_force(i.m*self.g)


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
		L = particles[0].pos-particles[1].pos #Find vector connecting the 2 particles
		#F= -k * (|L|- relaxed) * normal
		F = -particles[2] * (L.mag() - particles[3]) * L.norm()
		particles[0].add_force(F)
		particles[1].add_force(-F)


	def mouse_pull(self,particle):
		"""Force applied by selecting particle with mouse"""
		(pick_x,pick_y) = pygame.mouse.get_pos()
		mouse_pos = vect2d(pick_x,pick_y)
		self.selected.pos.x = mouse_pos.x
		self.selected.pos.y = mouse_pos.y
		"""dx = self.selected.pos - mouse_pos
		F = -self.mouse_force*dx - self.selected.vel*self.damping
		self.selected.add_force(F)"""



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
						self.selected.original_color = self.selected.color
						self.selected.color = (255,0,0)   #Change color
						self.induction(self.particle_list, picked)
			if event.type == pygame.MOUSEBUTTONUP and self.selected != None:
				self.selected.color = self.selected.original_color  #change color back to original
				self.force_que.remove(['mouse_pull',[self.selected]])   #Remove force from force_que
				self.selected = None    #No particles are selected

	def induction(self, particleList, position):

		"""Faraday's Law"""
		N = len(particleList[1].wire) #Number of Wires
		B = vect2d(position.x, position.y) #Magnetic Field Vector
		A = vect2d(position.x, position.y + particleList[0].height/2) #Parallel Vector to B
		PhiI = (B.x*A.x) + (B.y*A.y) #Phi's Equation (Initial) Needs to be DOT product NEEDS WORK
		PhiF = (B.x*A.x) + (B.y*A.y) #Phi's Equation (Final) Needs to be DOT product NEEDS WORK
		Phidx = PhiF - PhiI #Change in Phi
		E = -N*(Phidx/self.dt) #Number of wires*(Change in Phi / Change in time)

		print'E = {} N = {} B = {} A = {} PhiI = {} PhiF = {} Phidx = {}'.format(E, N, B, A, PhiI, PhiF, Phidx)
		
		"""Lenz's Law"""
		q = 1 #charge of particle
		EF = -1 #Electric Field
		v = B #original vector (Affected by Force)
		MF = B #magnetic Field Vector
		F = q*(EF + (v*MF)) #Force Equation

		print'F = {}'.format(F)


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
				rc = contact[0].getContactGeomParams()[0]
				n_hat = contact[0].getContactGeomParams()[1]
				dx = contact[0].getContactGeomParams()[2]
				#self.resolve_collision_ode(p[0], p[1], vect2d(rc[0], rc[1]), vect2d(n_hat[0], n_hat[1]), dx, contact)

	def bounding_sphere(self,particle1,particle2):
		"""Check for overlap using bounding spheres"""
		if ((particle2.pos - particle1.pos).mag()) < particle1.size + particle2.size:
			if particle1.shape =="circle" and particle2.shape == "circle":
				self.has_collided.append([particle1, particle2])
			elif particle1.shape == "rectangle" and particle2.shape == "rectangle":
				#Do nothing here, will be implemented at a  later date
				print("Rectangles have collided")
			elif particle1.shape == "rectangle" or particle2.shape == "rectangle":
				#Assignment 9 implementation
				#Find which particle is the rectangle
				rect = particle1
				part = particle2
				if particle2.shape == "rectangle":
					rect = particle2
					part = particle1

				#Find relative location of the particle to the rectangle
				r_rel = rect.pos - part.pos
				xb = r_rel * rect.local_x_axis.norm()
				yb = r_rel * rect.local_y_axis.norm()

				#Check if particle is within bounds of rectangle
				if xb+part.size >= -rect.width/2 and xb-part.size <= rect.width/2 and yb+part.size >= -rect.height/2 and yb - part.size <= rect.height/2:
					self.has_collided.append([particle1, particle2])


	def resolve_collision_ode(self, p1, p2, r_c, n_hat, dx, contact):

		#point particle vs point particle collision
		if p1.shape == "circle" and p2.shape == "circle":

			denom = p1.m_inv + p2.m_inv
			numer = (1 + (p1.cor * p2.cor)) * (p2.vel * n_hat - p1.vel * n_hat)

			dv1 = p1.m_inv * numer/denom
			dv2 = -p2.m_inv * numer/denom

			p1.vel += dv1 * n_hat
			p2.vel += dv2 * n_hat

			p1.pos += n_hat * dx * p1.m_inv/(p1.m_inv+p2.m_inv)
			p2.pos += n_hat * dx * -p2.m_inv/(p1.m_inv+p2.m_inv)

		#point particle vs rigid body collision
		elif (p1.shape == "circle" and p2.shape == "rectangle") or (p1.shape == "rectangle" and p2.shape == "circle"):
			#Find which particle is the rectangle
			rect = p1
			part = p2
			if p2.shape == "rectangle":
				rect = p2
				part = p1

			r_c -= rect.pos
			r_c_perp = vect2d(-r_c.y, r_c.x)
			#velocity of the contact point
			#v_c = vect2d(rect.vel.x + rect.omega*r_c_perp.x, rect.vel.y + rect.omega * r_c_perp.y)
			v_c = rect.vel + (rect.omega * r_c_perp)
			#Relative velocity between contact point and particle
			v_rel = v_c - part.vel
			v_close = v_rel * n_hat
			
			numer = (1 + (rect.cor * part.cor)) * v_close
			denom = rect.m_inv + part.m_inv
			I = rect.m/12 * (rect.width**2 + rect.height**2)
			J = numer/(denom + (((r_c_perp * n_hat)**2) / I))
				
			dv1 = (-J / rect.m) * n_hat
			dv2 = (J / part.m) * n_hat
			dw = ((J/I) * n_hat) * -r_c_perp

			rect.vel += dv1
			part.vel += dv2
			rect.omega += dw

			dp1 = n_hat * dx * p1.m_inv/(p1.m_inv + p2.m_inv)
			dp2 = n_hat * dx * -p2.m_inv/(p1.m_inv + p2.m_inv)

			p1.pos += n_hat * dx * p1.m_inv/(p1.m_inv + p2.m_inv)
			p2.pos += n_hat * dx * -p2.m_inv/(p1.m_inv + p2.m_inv)

		#rigid body vs rigid body collision
		elif p1.shape == "rectangle" and p2.shape == "rectangle":
			
			dv1 = vect2d(0,0)
			dv2 = vect2d(0,0)
			dw1 = 0
			dw2 = 0
			contact_num = int(len(contact)/2)#The number of contact points in this collision
			i = 0
			for j in range(contact_num):

				#Grab this particular contact point we're going to handle
				rc = contact[i].getContactGeomParams()[0]
				n_hat = contact[i].getContactGeomParams()[1]
				dx = contact[i].getContactGeomParams()[2]
				r_c = vect2d(rc[0], rc[1])
				n_hat = vect2d(n_hat[0], n_hat[1])

				#Get relative contact points for each rectangle along with the contact perps
				r_c1 = r_c - p1.pos
				r_c2 = r_c - p2.pos
				r_c_perp1 = vect2d(-r_c1.y, r_c1.x)
				r_c_perp2 = vect2d(-r_c2.y, r_c2.x)

				#Velocities of the contact points
				v_c1 = p1.vel + (p1.omega * r_c_perp1)
				v_c2 = p2.vel + (p2.omega * r_c_perp2)

				#Relative and closing velocities
				v_rel = v_c1 - v_c2
				v_close = v_rel * n_hat

				#Calculate the force
				numer = (1+(p1.cor * p2.cor)) * v_close
				denom = p1.m_inv + p2.m_inv
				I1 = p1.m/12 * (p1.width**2 + p1.height**2)
				I2 = p2.m/12 * (p2.width**2 + p2.height**2)
				J = numer/(denom + ((r_c_perp1 * n_hat)**2/I1) + ((r_c_perp2 * n_hat)**2/I2))
				
				#Sum the results for each each collision
				dv1 += (-J / p1.m) * n_hat
				dv2 += (J / p2.m) * n_hat
				dw1 += ((J/I1) * n_hat) * -r_c_perp1
				dw2 += ((J/I2) * n_hat) * r_c_perp2

				#Move foward in the list of contact points by 2
				i+=2

			#Take the average of all the data gathered
			dv1 /= contact_num
			dv2 /= contact_num
			dw1 /= contact_num
			dw2 /= contact_num
			#Add on new values
			p1.vel += dv1
			p1.omega += dw1
			p2.vel += dv2
			p2.omega += dw2

			#Resolve interpenetration
			p1.pos += n_hat * dx * p1.m_inv/(p1.m_inv + p2.m_inv)
			p2.pos += n_hat * dx * -p2.m_inv/(p1.m_inv + p2.m_inv)


	def resolve_collision(self):
		"""Changes velocities of point particles that have collided"""
		for i in self.has_collided:

			#Two particles have collided with eachother
			if i[0].shape == "circle" and i[1].shape == "circle":
				#Find parts of the equation for later use
				hat =  self.collision_normal(i[0], i[1])
				denom = i[0].m_inv + i[1].m_inv
				numer = (1 + (i[0].cor * i[1].cor)) * (i[1].vel * hat - i[0].vel * hat)

				#Calculate the change in velocities
				dv1 =  i[0].m_inv * numer / denom #* (i[0].m*i[1].m)
				dv2 = -i[1].m_inv * numer / denom #* (i[0].m*i[1].m)
				
				#Add change in velocity to current velocities
				i[0].vel += dv1 * hat
				i[1].vel += dv2 * hat

				self.resolve_interpenetration(i[0], i[1], hat)

			#Two rectangles have collided with eachother.  Currently does nothing
			elif i[0].shape =="rectangle" and i[1].shape == "retcangle":
				print("Recatngle collision resolution")

			#A rectangle has collided with another shape
			elif i[0].shape == "rectangle" or i[1].shape == "rectangle":
				#Find which particle is the rectangle
				rect = i[0]
				part = i[1]
				if i[1].shape == "rectangle":
					rect = i[1]
					part = i[0]

				#Find collision point
				#Find relative location of the particle to the rectangle
				r_c = vect2d(0,0)
				hat = vect2d(0,0)
				dx = 0
				r_rel_world = part.pos-rect.pos
				r_rel = vect2d(r_rel_world * rect.local_x_axis, r_rel_world * rect.local_y_axis)

				#Collided vertically
				if r_rel.x > -rect.width/2 and r_rel.x < rect.width/2:
					r_c.x = r_rel.x
					print "top/bot collision"

					if r_rel.y > rect.height/2:
						r_c.y = rect.height/2
						hat = rect.local_y_axis
						dx = part.size + rect.height/2 - r_rel.y
					else:
						r_c.y = -rect.height/2
						hat = -rect.local_y_axis
						dx = part.size + rect.height/2 + r_rel.y

				#We're contacting on the side of the rectangle
				elif r_rel.y > -rect.height/2 and r_rel.y < rect.height/2:
					print "side collision"
					r_c.y = r_rel.y

					if r_rel.x > rect.width/2:
						r_c.x = rect.width/2
						hat = rect.local_x_axis
						dx = part.size + rect.width/2 - r_rel.x
					else:
						r_c.x = -rect.width/2
						hat = -rect.local_x_axis
						dx = part.size + rect.width/2 + r_rel.x

				#Only option left is that we're colliding with a corner		
				else:
					print "corner collision"
					r_c = vect2d(rect.width/2, rect.height/2)

					if r_rel.x <= -rect.width/2:
						r_c.x *= -1
						flag_x = True
					if r_rel.y <= -rect.height/2:
						r_c.y *= -1
						flag_y = True

					#Check the flags to find the collision normal going to a corner
					hat = -(i[1].pos - i[0].pos).norm()
					dx = part.size + rect.size - (rect.pos - part.pos).mag()
					
				#Make sure collision normal is actually normalized
				hat = hat.norm()
					
				#contact point in world coordinates
				r_c = r_c.x * rect.local_x_axis.norm() + r_c.y * rect.local_y_axis.norm()
				r_c_perp = vect2d(-r_c.y, r_c.x)
				#velocity of the contact point
				#v_c = rect.vel + r_c_perp * rect.omega
				v_c = vect2d(rect.vel.x - rect.omega*r_c_perp.y, rect.vel.y + rect.omega * r_c_perp.x)
				#Relative velocity between contact point and particle
				v_rel = v_c - part.vel
				v_close = v_rel * hat

				numer = (1 + (rect.cor * part.cor)) * v_close
				denom = rect.m_inv + part.m_inv
				I = rect.m/12 * (rect.width**2 + rect.height**2)
				J = numer/(denom + ((r_c_perp * hat) / I))
				print "J: ", J
				
				dv1 = (-J / rect.m)
				dv2 = (J / part.m)
				dw = ((J * hat)/I) * -r_c_perp

				#print "dw: ", dw

				rect.vel += dv1 * hat
				part.vel += dv2 * hat
				rect.omega += dw * 0.01

				print "hat: ", hat
				print "dx: ", dx
				i[0].pos += hat * dx * i[0].m_inv/(i[0].m_inv + i[1].m_inv)
				i[1].pos += hat * dx * -i[1].m_inv/(i[0].m_inv + i[1].m_inv)

				print"-------------------------------------"


	def collision_normal(self,particle1,particle2):
		"""Returns the collision normal vector"""
		return -(particle1.pos-particle2.pos).norm()



	def resolve_interpenetration(self,particle1,particle2,hat):
		"""Move objects that overlap far enough apart that they don't continue to collide"""
		dx = particle1.size + particle2.size - (particle1.pos - particle2.pos).mag()

		particle1.pos += hat * dx * -particle1.m_inv/(particle1.m_inv+particle2.m_inv)
		particle2.pos += hat * dx * particle2.m_inv/(particle1.m_inv+particle2.m_inv)


	def boundary(self):
		"""Checks to see if part of a particle leaves the screen.  Should shift particle back on to screen and reverse component of velocity
		This method has been moved from the world class.  It makes more sense to have it be part of the collision engine
		"""
		for particle in self.particle_list:
			if particle.pos.x + particle.size >= self.world.width:
				particle.pos.x = self.world.width - particle.size
				particle.vel.x = -self.cor*particle.cor*particle.vel.x
			if particle.pos.x - particle.size <= 0:
				particle.pos.x = 0 + particle.size
				particle.vel.x = -self.cor*particle.cor*particle.vel.x
			if particle.pos.y + particle.size >= self.world.height:
				particle.pos.y = self.world.height - particle.size
				particle.vel.y = -self.cor*particle.cor*particle.vel.y
			if particle.pos.y - particle.size <= 0:
				particle.pos.y = 0 + particle.size
				particle.vel.y = -self.cor*particle.cor*particle.vel.y


	def check_rod_constraint(self,part0, part1,L):
		"""NEW: Check to see if particles are too close or too far apart.  Call fix_rod_constraint if Needed
		L is the length of the rod.  To get the length L you will need to use part0.rod_connection[part1].  This uses the rod_connection dictionary to look
		up the associated value, which is the rod length L.
		You will need to calculate the vector along the rod (hat) and
		the current distance between the particles (dist) to pass to fix_rod_constraint.
		"""
		dist = (part0.pos - part1.pos).mag()
		if dist > L or dist < L:
			self.fix_rod_constraint(part0, part1, L, self.collision_normal(part0, part1), dist) 

		

	def fix_rod_constraint(self,part0,part1,L,hat,dist):
		"""NEW: Move particle and change velocity so rod constraints are maintained
		L is length of rod, hat is the vector pointing along the rod, and dist is the
		current distance between the particles.  You can use the same impulse equations
		used during collisions to change the velocities.  To move the particles use the 
		interpentration resolution formula.
		"""
		#Calculate resulting velocities
		denom = part0.m_inv + part1.m_inv
		numer = (1 + (part0.cor * part1.cor)) * (part1.vel * hat - part0.vel * hat)

		#Calculate the change in velocities
		dv1 =  part0.m_inv * numer / denom #* (i[0].m*i[1].m)
		dv2 = -part1.m_inv * numer / denom #* (i[0].m*i[1].m)
			
		#Add change in velocity to current velocities
		part0.vel += dv1 * hat
		part1.vel += dv2 * hat


		#Resolve the interpenetration with the rod
		dx = L - dist

		part0.pos += hat * dx * -part0.m_inv/(part0.m_inv+part1.m_inv)
		part1.pos += hat * dx * part1.m_inv/(part0.m_inv+part1.m_inv)

		


if __name__ == '__main__': 
	#The following is just to test your code out and make sure it works
	earth = world()

	earth.set_numerical('euler')  #Change to the Velocity Verlet method

	part_list = []

	magnet = magnet(vect2d(350,200))
	part_list.append(magnet)
	earth.add_particle(magnet)

	coil = coil(vect2d(600, 200))

	earth.add_particle(coil)
	#earth.new_force('const_grav',part_list)   #Particles experience gravity
	#Leave the drag force in for use with moving the magnet around
	earth.new_force('drag',part_list)   #particles also experience drag
	
	"""
	for i in range(circle_num):
		for j in range(circle_num-1):
			earth.add_rod(part_list[i], part_list[j])
	"""
	###
	#While loop
	###

	while earth.running:   
		time_passed = earth.clock.tick(120)   #Add pygame.time.Clock() to world setup.  This will allow you to set fps
		earth.update()
