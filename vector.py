# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10  2015

Created by Todd Zimmerman for PHYS-360 Spring 2015 at UW-Stout

More advanced vect2d class, using + and - as well as scalar multiplication and division
Save assignment as program_homework_2.py and be sure to add your name to the author_name variable below
"""

from __future__ import division #imports non-integer division.  Without this you get 1/2 = 0, not 1/2 = 0.5
from math import sqrt  #Imports the square root function in case you need it


######
author_name = "Aaron Aumann"  #Include your name as a string to get credit for completing assignment
######

class vect2d():
    def __init__(self,x=1,y=0):
        """2D Vector class (x,y)
        
        x and y variables have default values so that if you don't
        specify values for those variables the default values will be used
        """        
        self.x = float(x)   #Components are set as floating point numbers because
        self.y = float(y)   # Python uses integer division by default where 5/4 = 1 (i.e. decimal dropped in integer div)

    def __len__(self):
        return 2

        
    def __rep__(self):
        """Method used internally to call vectors, must return a string that looks like vect3d object
        
        """
        return "vect2d(%s,%s)" % (self.x,self.y)
        
    
    def __str__(self):
        """Method used by print command, returns a string
        
        Note about docstrings: If were to input the text after the greater-than signs
        you should get the output below those lines
        if your method is correct
        >>>A = vect2d(3,2))
        >>>print A
        (3,2)
        
        Note that %s is replaced by value of self.x or self.y
        """
        return "(%s, %s)" % (self.x, self.y)  


    def __eq__(self,other):
        """Checks to see if two vectors are equal to each other
        >>>A = vect2d(2,3)
        >>>B = vect2d(2,3)
        >>>print A == B
        True
        """
        if self.x == other.x and self.y == other.y:
            return True
        else:
            return False

    def __add__(self,other):
        """Add together vectors:  self + other"""
        #Add your code here
        return vect2d(self.x + other.x, self.y + other.y)

    __radd__ = __add__
    """__radd__ stands for right addition.  not typically needed"""

    def __iadd__(self,other):
        """In place addition: self += other is same as self = self + other"""
        #Add your code here
        return self + other

    def __sub__(self,other):
        """Subtract other from self: self - other"""
        #Add your code here
        return vect2d(self.x - other.x, self.y - other.y)

    def __rsub__(self,other):
        """Subtract self from other: other - self"""
        #Add your code here
        return vect2d(other.x - self.x, other.y - self.y)

    def __isub__(self,other):
        """In place subtration: self -= other is same as self = self - other"""
        #Add your code here
        return self - other

    def __mul__(self,other):
        """Multiply a vector (self) by a scalar or vector (other)

        Make sure other isn't a vect2d by using the following test:
            isinstance(other, vect2d)
        This will return True if other is an instance of vect2d.
        Have method return dot product if both are vectors.
        """
        if isinstance(other,vect2d) and isinstance(self, vect2d):
            return self.x*other.x + self.y*other.y
        else:
            return vect2d(other*self.x, other*self.y)

    def __rmul__(self,other):
        """Multiply a vector (self) by a scalar or vector (other) when vector is right-hand operand

        Must define right multiplication because this will
        be used if the scalar comes first and the vector second.
        Have method return dot product if both are vectors"""
        if isinstance(self,vect2d) and isinstance(other,vect2d):
            return self.x*other.x + self.y*other.y
        else:
            return vect2d(self.x*other, self.y*other)


    def _imul___(self,other):
        """In place multiplication: self *= other is same as self = self*other

        other must be a scalar for this method to work. 
        Return None if other is a vector.
        """
        #Add your code here
        return self * other


    def __div__(self,other):
        """Division of self by other: self/other

        Return None if other is a vector as well since
        division by a vector is not defined.
        """
        #Add your code here
        if isinstance(other, vect2d):
            return None
        else:#Other is a scalar
            if other != 0:#Avoid division by 0
                return vect2d(self.x / other, self.y / other)
            else:#If dividing by zero, return None
                return None

    def __rdiv__(self,other):
        """rdiv is only called if the right hand operand is a vector so must return None"""
        #Add your code here
        return None

    def __idiv__(self,other):
        """In place division: self/=other is same as self = self/other

        If other is vect2d then return None.
        """
        #Add your code here
        return self / other


    def __truediv__(self,other):
        """Division of self by other when using division from __future__: self/other

        True division is what Python calls whenever you import the
        division method from the __future__ module.  Code should be same
        as for previous division methods.
        Return None if other is a vector as well since
        division by a vector is not defined.
        """
        #Add your code here
        if isinstance(other, vect2d):
            return None
        else:#Other is a scalar
            if other != 0:#Avoid division by 0
                return vect2d(self.x / other, self.y / other)
            else:#If dividing by zero, return None
                return None

    def __rtruediv__(self,other):
        """rdiv is only called if the right hand operand is a vector so must return None"""
        #Add your code here
        return None

    def __itruediv__(self,other):
        """In place division: self/=other is same as self = self/other

        If other is vect2d then return None.
        """
        #Add your code here
        return self/other

    def __neg__(self):
        """Returns the negative of the vector: -A"""
        #Add your code here
        return self * -1
        
    def mag(self):
        """Find the magnitude of the vector
        
        >>>A = vect2d(3,4)
        >>>D = A.mag()
        >>>print D
        5
        """
        #Reuse code from previous assignment
        return sqrt(self.x*self.x + self.y * self.y)
        
    def mag2(self):
        """Find the magnitude squared

        >>>A = vect2d(3,4)
        >>>D = A.mag2()
        25
        """
        #Reuse code from previous assignment
        return self.x * self.x + self.y * self.y

        
    def norm(self):
        """Calculate the unit vector 
        
        The unit vector is called a normalized vector.  Normalizing
        is rescaling a vector so that its length (magnitude) is one
        >>>A = vect2d(2,8)
        >>>E = A.norm()
        >>>print E
        (0.2425, 0.9701)
        """
        #Reuse code from previous assignment
        #Avoid divsion by 0
        if self.mag() == 0:
            return vect2d(0,0)
            
        x = self.x / self.mag()
        y = self.y / self.mag()
        return vect2d(x, y)

    def dot(self,other):
        """Calculates the dot product of two vectors"""
        if isinstance(other,vect2d):
            return self.x*other.x + self.y*other.y
        else:
            return None
            
    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        else:
            raise IndexError("Index out of range "+str(key)+" for vec2d")


    def __setitem__(self, index, value):
        if index == 0:
            self.x = value
        elif index == 1:
            self.y = value
        else:
            raise IndexError("Index out of range "+str(key)+" for vec2d")

    def integ(self):
        """Return a vect2d with integer values

        This is useful for passing vectors as points in the pygame world
        """
        return vect2d(int(self.x), int(self.y))
           
          
"""
Vector 3d class
"""  
class vect3d():
    def __init__(self, x=1, y=0, z=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def length():
        return 2

    def __rep__(self):
        """Method used internally to call vectors, must return a string that looks like vect3d object
        
        """
        return "vect3d(%s,%s,%s)" % (self.x,self.y, self.z)
        
    
    def __str__(self):
        """Method used by print command, returns a string
        
        Note about docstrings: If were to input the text after the greater-than signs
        you should get the output below those lines
        if your method is correct
        >>>A = vect2d(3,2))
        >>>print A
        (3,2)
        
        Note that %s is replaced by value of self.x or self.y
        """
        return "(%s, %s,%s)" % (self.x, self.y,self.z)  


    def __eq__(self,other):
        """Checks to see if two vectors are equal to each other
        >>>A = vect2d(2,3)
        >>>B = vect2d(2,3)
        >>>print A == B
        True
        """
        if self.x == other.x and self.y == other.y and self.z == other.z:
            return True
        else:
            return False

    def __add__(self,other):
        """Add together vectors:  self + other"""
        #Add your code here
        return vect3d(self.x + other.x, self.y + other.y, self.z + other.z)

    __radd__ = __add__
    """__radd__ stands for right addition.  not typically needed"""

    def __iadd__(self,other):
        """In place addition: self += other is same as self = self + other"""
        #Add your code here
        return self + other

    def __sub__(self,other):
        """Subtract other from self: self - other"""
        #Add your code here
        return vect3d(self.x - other.x, self.y - other.y, self.z - other.z)

    def __rsub__(self,other):
        """Subtract self from other: other - self"""
        #Add your code here
        return vect3d(other.x - self.x, other.y - self.y, other.z - self.z)

    def __isub__(self,other):
        """In place subtration: self -= other is same as self = self - other"""
        #Add your code here
        return self - other

    def __mul__(self,other):
        """Multiply a vector (self) by a scalar or vector (other)

        Make sure other isn't a vect2d by using the following test:
            isinstance(other, vect2d)
        This will return True if other is an instance of vect2d.
        Have method return dot product if both are vectors.
        """
        if isinstance(other,vect3d) and isinstance(self, vect3d):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            return vect3d(other*self.x, other*self.y, other*self.z)

    def __rmul__(self,other):
        """Multiply a vector (self) by a scalar or vector (other) when vector is right-hand operand

        Must define right multiplication because this will
        be used if the scalar comes first and the vector second.
        Have method return dot product if both are vectors"""
        if isinstance(self,vect3d) and isinstance(other,vect3d):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            return vect3d(self.x*other, self.y*other, self.z*other)


    def _imul___(self,other):
        """In place multiplication: self *= other is same as self = self*other

        other must be a scalar for this method to work. 
        Return None if other is a vector.
        """
        #Add your code here
        return self * other


    def __div__(self,other):
        """Division of self by other: self/other

        Return None if other is a vector as well since
        division by a vector is not defined.
        """
        #Add your code here
        if isinstance(other, vect3d):
            return None
        else:#Other is a scalar
            if other != 0:#Avoid division by 0
                return vect3d(self.x / other, self.y / other, self.z / other)
            else:#If dividing by zero, return None
                return None

    def __rdiv__(self,other):
        """rdiv is only called if the right hand operand is a vector so must return None"""
        #Add your code here
        return None

    def __idiv__(self,other):
        """In place division: self/=other is same as self = self/other

        If other is vect2d then return None.
        """
        #Add your code here
        return self / other


    def __truediv__(self,other):
        """Division of self by other when using division from __future__: self/other

        True division is what Python calls whenever you import the
        division method from the __future__ module.  Code should be same
        as for previous division methods.
        Return None if other is a vector as well since
        division by a vector is not defined.
        """
        #Add your code here
        if isinstance(other, vect3d):
            return None
        else:#Other is a scalar
            if other != 0:#Avoid division by 0
                return vect3d(self.x / other, self.y / other, self.z / other)
            else:#If dividing by zero, return None
                return None

    def __rtruediv__(self,other):
        """rdiv is only called if the right hand operand is a vector so must return None"""
        #Add your code here
        return None

    def __itruediv__(self,other):
        """In place division: self/=other is same as self = self/other

        If other is vect2d then return None.
        """
        #Add your code here
        return self/other

    def __neg__(self):
        """Returns the negative of the vector: -A"""
        #Add your code here
        return self * -1
        
    def mag(self):
        """Find the magnitude of the vector
        
        >>>A = vect2d(3,4)
        >>>D = A.mag()
        >>>print D
        5
        """
        #Reuse code from previous assignment
        return sqrt(self.x*self.x + self.y * self.y + self.z * self.z)
        
    def mag2(self):
        """Find the magnitude squared

        >>>A = vect2d(3,4)
        >>>D = A.mag2()
        25
        """
        #Reuse code from previous assignment
        return self.x * self.x + self.y * self.y + self.z * self.z

        
    def norm(self):
        """Calculate the unit vector 
        
        The unit vector is called a normalized vector.  Normalizing
        is rescaling a vector so that its length (magnitude) is one
        >>>A = vect2d(2,8)
        >>>E = A.norm()
        >>>print E
        (0.2425, 0.9701)
        """
        #Reuse code from previous assignment
        #Avoid divsion by 0
        if self.mag() == 0:
            return vect3d(0,0)
            
        x = self.x / self.mag()
        y = self.y / self.mag()
        z = self.z / self.mag()
        return vect3d(x, y, z)

    def dot(self,other):
        """Calculates the dot product of two vectors"""
        if isinstance(other,vect2d):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            return None
            
    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError("Index out of range "+str(key)+" for vec3d")


    def __setitem__(self, index, value):
        if index == 0:
            self.x = value
        elif index == 1:
            self.y = value
        elif index == 2:
            self.z = value
        else:
            raise IndexError("Index out of range "+str(key)+" for vec2d")

    def integ(self):
        """Return a vect2d with integer values

        This is useful for passing vectors as points in the pygame world
        """
        return vect3d(int(self.x), int(self.y), int(self.z))

########################
#You don't need any of the following code, it is just to try out your methods

if __name__ == '__main__' :  #All of the following code only runs if you run this file as your main program
        
    A=vect2d(2,8)
    B=vect2d(-1,2)
    print "Number of elements in A is", len(A)
    print "Vector A equals ", A
    print "Vector B equals ", B
    print "A + B = ", A+B
    print "A - B = ", A-B
    print "3 * A = ", 3*A
    print "A * 3 = ", A*3
    print "A/3 =", A/3
    print "3/A =", 3/A
    print A.mag()
    print A.norm()
    print A.mag2()

    X = vect3d(2,8,10)
    Y = vect3d(-1,2,5)
    print "Vector X equals ", X
    print "Vector Y equals ", Y
    print "X + Y = ", X+Y
    print "X - Y = ", X-Y
    print "3 * X = ", 3*X
    print "X * 3 = ", X*3
    print "X/3 =", X/3
    print "3/X =", 3/X
    print X.mag()
    print X.norm()
    print X.mag2()
