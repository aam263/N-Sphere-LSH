import numpy as np
from scipy.special import gamma
import scipy.integrate as integrate
from scipy import optimize
from collections import deque
import itertools

import time
import random as rand
import matplotlib.pyplot as plt

#################################################################################################################################

class spherical_tree_node:

	def __init__(self,median):
		self.median = median
		self.bounds = None
		self.left_child = None
		self.right_child = None

class spherical_hash_tree:

	"replace print statements with something better?"
	"in the memory regions where this structure is viable, hash_store is slightly faster"
	"this structure only allows for dimensons of up to around 22 before it crashes because of memory"

	def __init__(self,dimensions):
		if dimensions > 0:
			self.dimensions = dimensions
			self.root = spherical_tree_node(np.pi)
			self.__initialize(dimensions,self.root)
			self.curr_id = 0
			self.bounds_id = dict()
		else:
			print("error, expand later")


	def __initialize(self,dimensions,root):
		if dimensions > 0:
			left_child = spherical_tree_node(np.pi/2)
			root.left_child = left_child
			self.__initialize(dimensions-1,left_child)
			right_child = spherical_tree_node(np.pi/2)
			root.right_child = right_child
			self.__initialize(dimensions-1,right_child)
		else:
			root.bounds = []

	
	def insert(self, bounds_vec):
		self.__insert_helper(bounds_vec,0,self.root)


	def __insert_helper(self, bounds_vec, curr_dim, curr_node):
		if curr_dim > -1 and curr_dim < self.dimensions:
			if len(bounds_vec[curr_dim]) == 2:
				if curr_node.left_child != None and bounds_vec[curr_dim][0] < curr_node.median:
					#print("left")
					self.__insert_helper(bounds_vec,curr_dim+1,curr_node.left_child)
				elif curr_node.right_child != None and bounds_vec[curr_dim][0] >= curr_node.median:
					#print("right")
					self.__insert_helper(bounds_vec,curr_dim+1,curr_node.right_child)
				else:
					print("child nodes do not exist.")
			else:
				print("curr_dim is not of len 2 it is invalid")
		else:
			if curr_node.bounds != None:
				curr_node.bounds.append(bounds_vec)
			else:
				print("failure in construction of bounds in leaf nodes")
			

	def get_hash(self, unit_vec):
		"test"
		#assume hash has been created
		#assume unit vec in spherical coordinates
		"what to do is string_id is None?"
		string_id = self.__get_hash_helper(0,self.root,unit_vec)
		if string_id == None: print("no hash found; none value returned")
		return string_id

		

	def __get_hash_helper(self, dimensions, curr_node, unit_vec):
		"test - seems to work for 1 dimension sphere"
		#returns none if bounds not found
		if dimensions > -1 and dimensions < self.dimensions:
			if unit_vec[dimensions] <= curr_node.median:
				string_id = self.__get_hash_helper(dimensions+1,curr_node.left_child,unit_vec)
			elif unit_vec[dimensions] > curr_node.median: 
				string_id = self.__get_hash_helper(dimensions+1,curr_node.right_child,unit_vec)
			else:
				print("error, bounds not of appropriate array dimensions")
		elif dimensions == self.dimensions: 
			for bounds_vec in curr_node.bounds:
				#returns the first possible hash by default (only relevant for things that fall exactly on the edges of the spherical hash)
				_id = self.__get_id(bounds_vec,unit_vec)
				if _id != None: 
					string_id = _id
					break
			if string_id == None: print("no hash found; None value returned")
		else:
			print("dimensions out of bounds")
		return string_id

	def __get_id(self, bounds_vec, unit_vec):
		#assigns a hash to a vector based on bounds, with hashes generated for different bounds as they appear
		"test"
		found = True
		string_id = ""
		for idx, bounds in enumerate(bounds_vec):
			#if its on the border put it in the first hash it appears in
			if len(bounds) == 2:
				if unit_vec[idx] >= bounds[0] and unit_vec[idx] <= bounds[1]: 
					if bounds in self.bounds_id:
						string_id += str(self.bounds_id[bounds]) + "-" if idx != len(bounds_vec)-1 else str(self.bounds_id[bounds])
					else:
						string_id += str(self.curr_id) + "-" if idx != len(bounds_vec)-1 else str(self.curr_id)
						self.bounds_id[bounds] = self.curr_id
						self.curr_id += 1
				else:
					found = False
					break
			else:
				print("bounds in hash_tree are invalid. error")
		return string_id if found == True else None

##########################################################################################################################################
##########################################################################################################################################


class hash_store:

	"continue testing the spherical coordinate transformation"
	"remove redundancy in generate spherical hash code"

	curr_id = 0

	def __init__(self, left_bound=None, right_bound=None):
		self.left_bound = None
		self.right_bound = None
		self.children = []
		self.id = hash_store.curr_id
		hash_store.curr_id += 1

		if left_bound is not None:
			self.left_bound = left_bound
		if right_bound is not None:
			self.right_bound = right_bound

	def add(self,child):
		self.children.append(child)

	def get_children(self):
		return self.children

	def get_id(self):
		return str(self.id)

	
	@staticmethod
	def print_hash_bounds(root):
		hash_store.__print_hash_bounds_helper(root,"","")

	@staticmethod
	def __print_hash_bounds_helper(root, bound_string, id_string):
		if len(root.get_children()) > 0:
			for child in root.get_children():
				temp1 = "[{},{}]".format(child.left_bound,child.right_bound)+bound_string
				temp2 = child.get_id() + "-" + id_string 
				hash_store.__print_hash_bounds_helper(child,temp1,temp2)
		else:
			#if root has at least 1 child then the below is always safe
			if len(id_string) > 1:
				print(bound_string, id_string[0:-1])
			else:
				print("hash root has no children, and therefore no bounds or ids.")


	@staticmethod
	def get_hash(vector, root):
		#root must be the root defined by not having any bounds
		"not working in 1 dimension"
		hash_list = deque()
		i = 0
		hash_store.__get_hash_helper(vector,root,i,hash_list)
		hash_id = ""
		flag = False
		for idv_hash in hash_list:
			hash_id += (idv_hash + "-" if flag is False else idv_hash)
			flag = True
		return hash_id

	@staticmethod
	def __get_hash_helper(vector,root,i,hash_list):
		#assert that vector has the same dimesnions as the depth of the tree
		#if vector sits on borders of different regions then it falls in the hash that it first appears in
		if i < len(vector):
			for child in root.get_children():
				if vector[-i-1] >= child.left_bound and vector[-i-1] <= child.right_bound:
					#hash_list.append_left[child.left_bound,child.right_bound])
					hash_list.appendleft(child.get_id())
					hash_store.__get_hash_helper(vector,child,i+1,hash_list)



###########################################################################################################################################

class spherical_hash:


	@staticmethod
	def nsphere_area(dimensions):
		#d is number of dimensions
		if dimensions >= 0 and isinstance(dimensions,int):
			return 2*(np.pi**((dimensions+1.0)/2))/gamma(((dimensions+1.0)/2))
		else:
			print("Expression only valid for nonnegative integer dimensions; dimension={} is not supported".format(dimensions))


	@staticmethod
	def spherical_cap_area(dimensions, theta):
		#dimensions is dimesnions of the unit sphere, not cartesian dimensions
		#theta is radial angle of spherical cap
		#returns area of a spherical cap of radial angle theta
		#only ever used for dimesions > 1
		if dimensions > 0 and isinstance(dimensions,int) and theta >= 0 and theta <= np.pi:
			integrand = lambda x: np.sin(x)**(dimensions-1)
			return spherical_hash.nsphere_area(dimensions-1)*(integrate.quad(integrand, 0, theta))[0]
		elif dimensions <= 0 or not isinstance(dimensions,int):
			print("Expression only valid for positive integer dimensions; dimensions={} is not supported.".format(dimensions))
		elif theta < 0 or theta > np.pi:
			print("Expression defined for theta in [0,Pi]; theta={} is not in the specified range.".format(theta))


	@staticmethod
	def spherical_cap_angle(dimensions, area):
		#dimensions is number of dimensions of unit sphere, not cartesian dimensions
		#area is the area of this spherical cap
		#returns the angle theta s.t. the spherical cap of radial angle theta has the specified area
		area_bound = spherical_hash.nsphere_area(dimensions)
		if dimensions > 0 and isinstance(dimensions,int) and area >= 0 and area <= area_bound:
			function = lambda x: (spherical_hash.spherical_cap_area(dimensions, x)-area)
			return optimize.brentq(function, 0, np.pi)
		elif dimensions <= 0 or not isinstance(dimensions,int):
			print("Expression only valid for positive integer dimensions; dimensions={} is not supported.".format(dimensions))
		elif area < 0 or area > area_bound:
			print("Expression only valid for areas in [0,{}]; area={} is not in the specified range.".format(area_bound,area))


	@staticmethod
	def generate_spherical_hash_1(dimensions, partitions):
		#dimensions is the dimesnions of the unit sphere (cartesian dimensions is 1+dimensions)
		#partitions is how many equal area segments you want
		root = hash_store()
		spherical_hash.__generate_spherical_hash_helper_1(dimensions, partitions, root)
		return root


	@staticmethod
	def __generate_spherical_hash_helper_1(dimensions, partitions, root):
		#should return an ndlist or array
		"assert that root is instance of hashstore"
		if dimensions > 0 and isinstance(dimensions,int) and partitions > 0:
			desired_area = spherical_hash.nsphere_area(dimensions)/partitions
			if dimensions == 1:
				#is a partitioning of a circle
				for i in range(1,partitions+1):
					left_bound = (i-1)*(2*np.pi/partitions)
					right_bound = i*(2*np.pi/partitions)
					curr_root = hash_store(left_bound=left_bound,right_bound=right_bound)
					root.add(curr_root)
				
			elif dimensions != 1 and partitions == 1:
				#cap angle is pi and the single region is the entire sphere
				curr_root = hash_store(left_bound=0, right_bound=np.pi)
				spherical_hash.__generate_spherical_hash_helper_1(dimensions-1,partitions,curr_root)
				root.add(curr_root)

			else:
				cap_angle = spherical_hash.spherical_cap_angle(dimensions,desired_area)
				if partitions == 2:
					#just the north and the south polar caps
					upper_half = hash_store(left_bound=0,right_bound=cap_angle)
					spherical_hash.__generate_spherical_hash_helper_1(dimensions-1,partitions-1,upper_half)
					root.add(upper_half)
					lower_half = hash_store(left_bound=cap_angle,right_bound=np.pi)
					spherical_hash.__generate_spherical_hash_helper_1(dimensions-1,partitions-1,lower_half)
					root.add(lower_half)
				else:
					#if there is at least one collar
					upper_cap_root = hash_store(left_bound=0,right_bound=cap_angle)
					spherical_hash.__generate_spherical_hash_helper_1(dimensions-1,1,upper_cap_root)
					root.add(upper_cap_root)
					#########################################################################
					#deal with collars
					ideal_num_collars = max(1,int(np.floor(((np.pi-2*cap_angle)/(desired_area**(1.0/dimensions)))+0.5)))
					ideal_fitting_collar_angle = (np.pi-2*cap_angle)/ideal_num_collars
					curr_err = 0 # running sum of the error between ideal number of zones per collar and real number of zones per collar
					curr_area = desired_area
					real_collar_colatitude_1 = cap_angle
					for collar_num in range(1,ideal_num_collars+1):
						#collar is identified by the region inbetween colatitude_2 and colatitude_1
						ideal_collar_colatitude_1 = cap_angle + (collar_num-1)*ideal_fitting_collar_angle	
						ideal_collar_colatitude_2 = ideal_collar_colatitude_1 + ideal_fitting_collar_angle	
						num_ideal_collar_zones = (spherical_hash.spherical_cap_area(dimensions,ideal_collar_colatitude_2)-spherical_hash.spherical_cap_area(dimensions,ideal_collar_colatitude_1))/desired_area
						num_real_collar_zones = int(np.floor(num_ideal_collar_zones+curr_err+0.5))
						curr_err += (num_ideal_collar_zones-num_real_collar_zones)
						curr_area += num_real_collar_zones*desired_area
						real_collar_colatitude_2 = spherical_hash.spherical_cap_angle(dimensions,curr_area)
						curr_root = hash_store(left_bound=real_collar_colatitude_1,right_bound=real_collar_colatitude_2)
						spherical_hash.__generate_spherical_hash_helper_1(dimensions-1,num_real_collar_zones,curr_root)
						root.add(curr_root)
						real_collar_colatitude_1 = real_collar_colatitude_2
					#########################################################################
					lower_cap_root = hash_store(left_bound=real_collar_colatitude_2,right_bound=np.pi)
					spherical_hash.__generate_spherical_hash_helper_1(dimensions-1,1,lower_cap_root)
					root.add(lower_cap_root)

		elif dimensions	<= 0 or not isinstance(dimensions,int):
			print("Expression only valid for positive integer dimensions; dimensions={} is not supported.".format(dimensions))
		elif partitions	<= 0:
			print("Expression only valid if partitions is a positive integer; partitions={} is not supported".format(partitions))
		
##############################################################################################################################################################
	@staticmethod
	def generate_spherical_hash_2(dimensions, partitions):
		#dimensions is the dimesnions of the unit sphere (cartesian dimensions is 1+dimensions)
		#partitions is how many equal area segments you want
		curr_bounds = deque()
		hash_tree = spherical_hash_tree(dimensions)
		spherical_hash.__generate_spherical_hash_helper_2(dimensions, partitions, curr_bounds, hash_tree)
		return hash_tree


	@staticmethod
	def __generate_spherical_hash_helper_2(dimensions, partitions, curr_bounds, hash_tree):
		#should return an ndlist or array
		"assert that hash_tree is instance of spherical_hash_tree"
		if dimensions > 0 and isinstance(dimensions,int) and partitions > 0:
			desired_area = spherical_hash.nsphere_area(dimensions)/partitions
			if dimensions == 1:
				#is a partitioning of a circle
				for i in range(1,partitions+1):
					left_bound = (i-1)*(2*np.pi/partitions)
					right_bound = i*(2*np.pi/partitions)
					curr_bounds.appendleft((left_bound,right_bound))###########111
					hash_tree.insert(list(curr_bounds)) 
					curr_bounds.popleft()
				
			elif dimensions != 1 and partitions == 1:
				#cap angle is pi and the single region is the entire sphere
				curr_bounds.appendleft((0,np.pi)) ##############111
				spherical_hash.__generate_spherical_hash_helper_2(dimensions-1,partitions,curr_bounds,hash_tree)
				curr_bounds.popleft() # not necessary

			else:
				cap_angle = spherical_hash.spherical_cap_angle(dimensions,desired_area)
				if partitions == 2:
					#just the north and the south polar caps
					curr_bounds.appendleft((0,cap_angle))############111
					spherical_hash.__generate_spherical_hash_helper_2(dimensions-1,partitions-1,curr_bounds,hash_tree)
					curr_bounds.popleft()
					curr_bounds.appendleft((cap_angle,np.pi))#############111
					spherical_hash.__generate_spherical_hash_helper_2(dimensions-1,partitions-1,curr_bounds,hash_tree)
					curr_bounds.popleft()
				else:
					#if there is at least one collar
					curr_bounds.appendleft((0,cap_angle))#################111
					spherical_hash.__generate_spherical_hash_helper_2(dimensions-1,1,curr_bounds,hash_tree)
					curr_bounds.popleft()
					#########################################################################
					#deal with collars
					ideal_num_collars = max(1,int(np.floor(((np.pi-2*cap_angle)/(desired_area**(1.0/dimensions)))+0.5)))
					ideal_fitting_collar_angle = (np.pi-2*cap_angle)/ideal_num_collars
					curr_err = 0 # running sum of the error between ideal number of zones per collar and real number of zones per collar
					curr_area = desired_area
					real_collar_colatitude_1 = cap_angle
					for collar_num in range(1,ideal_num_collars+1):
						#collar is identified by the region inbetween colatitude_2 and colatitude_1
						ideal_collar_colatitude_1 = cap_angle + (collar_num-1)*ideal_fitting_collar_angle	
						ideal_collar_colatitude_2 = ideal_collar_colatitude_1 + ideal_fitting_collar_angle	
						num_ideal_collar_zones = (spherical_hash.spherical_cap_area(dimensions,ideal_collar_colatitude_2)-spherical_hash.spherical_cap_area(dimensions,ideal_collar_colatitude_1))/desired_area
						num_real_collar_zones = int(np.floor(num_ideal_collar_zones+curr_err+0.5))
						curr_err += (num_ideal_collar_zones-num_real_collar_zones)
						curr_area += num_real_collar_zones*desired_area
						real_collar_colatitude_2 = spherical_hash.spherical_cap_angle(dimensions,curr_area)
						curr_bounds.appendleft((real_collar_colatitude_1,real_collar_colatitude_2))###########111
						spherical_hash.__generate_spherical_hash_helper_2(dimensions-1,num_real_collar_zones,curr_bounds,hash_tree)
						curr_bounds.popleft()
						real_collar_colatitude_1 = real_collar_colatitude_2
					#########################################################################
					curr_bounds.appendleft((real_collar_colatitude_2,np.pi))#############111
					spherical_hash.__generate_spherical_hash_helper_2(dimensions-1,1,curr_bounds,hash_tree)
					curr_bounds.popleft() # not necessary

		elif dimensions	<= 0 or not isinstance(dimensions,int):
			print("Expression only valid for positive integer dimensions; dimensions={} is not supported.".format(dimensions))
		elif partitions	<= 0:
			print("Expression only valid if partitions is a positive integer; partitions={} is not supported".format(partitions))
		

#################################################################3#################################################################################################
	
	@staticmethod
	def to_spherical_coordinates(vector):
		"seems to work, more intensive testing probably a good idea"
		#works for the map (psi: [0,2*pi]x[0,pi]^{d-1} -> s^{d} contained in R^{d+1})
		norm = np.linalg.norm(vector)
		normed_vec = [x/norm for x in vector]
		coordinates = deque()
		cartesian_dim = len(vector)
		if cartesian_dim > 1:
			curr_divisor = 1
			while cartesian_dim > 1:
				if curr_divisor != 0:
					curr_val = normed_vec[cartesian_dim-1]/curr_divisor
					if cartesian_dim == 2:
						theta = np.arcsin(curr_val)
						if normed_vec[cartesian_dim-2] < 0:
							theta = np.pi - theta
						elif normed_vec[cartesian_dim-2] > 0 and normed_vec[cartesian_dim-1] < 0:
							theta = 2*np.pi + theta
					else:
						theta = np.arccos(curr_val) 
				else:
					#when a pole is encountered we can choose arbitrary value is appropriate range for angles. Choose 0.
					theta = 0.0
					print("Pole encountered: arbitrary values set to 0.")
				coordinates.appendleft(theta) 
				curr_divisor *= (np.sin(theta) if cartesian_dim != 2 else 1)
				cartesian_dim -= 1
			return coordinates
		else:
			print("Expression only valid for vectors with dimensions>1; dimensions={} is not valid.".format(len(vector)))

################################################################################################
	def __init__(self,dimensions,partitions):
		#dimesnions is the dimensions of the unit sphere 
		self.__root = spherical_hash.generate_spherical_hash_1(dimensions,partitions)
		#self.hash_tree = spherical_hash.generate_spherical_hash_2(dimensions,partitions)

	def get_hash_1(self, vector, left_bound, right_bound):
		#get the hash for the current vectors in the dimensions in range(left_bound,right_bound)
		coordinates = spherical_hash.to_spherical_coordinates(vector[left_bound:right_bound])
		return hash_store.get_hash(vector,self.__root)

#	def get_hash_2(self, vector, left_bound, right_bound):
#		#get the hash for the current vectors in the dimensions in range(left_bound,right_bound)
#		coordinates = spherical_hash.to_spherical_coordinates(vector[left_bound:right_bound])
#		return self.hash_tree.get_hash(coordinates)




##########################################################################################################################################################
	
	

				




if __name__ == "__main__":
	"need to fix get_hash_1 for 1 dimension"
	dimensions = [dim for dim in range(1,250)]
	partitions1 = 10
	partitions2 = 50
	partitions3 = 100
	partitions4 = 150
	times1 = []
	times2 = []
	times3 = []
	times4 = []
	
	for dim in dimensions:
		#spherical dim
		
		h = spherical_hash(dim,partitions1)
		print(dim)
		t1 = time.time()
		for i in range(0,1000):
			vec = np.random.uniform(low=0, high=5, size=(dim+1,))
			x = h.get_hash_1(vec,0,len(vec))
		times1.append((time.time()-t1)/1000)

	for dim in dimensions:
		#spherical dim
		
		h = spherical_hash(dim,partitions2)
		print(dim)
		t2 = time.time()
		for i in range(0,1000):
			vec = np.random.uniform(low=0, high=5, size=(dim+1,))
			x = h.get_hash_1(vec,0,len(vec))
		times2.append((time.time()-t2)/1000)

	for dim in dimensions:
		#spherical dim
		
		h = spherical_hash(dim,partitions3)
		print(dim)
		t1 = time.time()
		for i in range(0,1000):
			vec = np.random.uniform(low=0, high=5, size=(dim+1,))
			x = h.get_hash_1(vec,0,len(vec))
		times3.append((time.time()-t1)/1000)

	for dim in dimensions:
		#spherical dim
		
		h = spherical_hash(dim,partitions4)
		print(dim)
		t1 = time.time()
		for i in range(0,1000):
			vec = np.random.uniform(low=0, high=5, size=(dim+1,))
			x = h.get_hash_1(vec,0,len(vec))
		times4.append((time.time()-t1)/1000)
		
		#t2 = time.time()
		#for i in range(0,100):
		#	vec = np.random.uniform(low=0, high=5, size=(dim+1,))
		#	x = h.get_hash_2(vec,0,len(vec))
		#times2.append((time.time()-t2)/100)

		del h

	#ax1.scatter(dimensions, times2, s=10, c='r', marker="o", label='second')
	plt.plot(dimensions,times1,"r",label="p = 25")
	plt.plot(dimensions,times2,"g",label="p = 50")
	plt.plot(dimensions,times3,"b",label="p = 100")
	plt.plot(dimensions,times4,"y",label="p = 200")
	plt.xlabel("dim")
	plt.ylabel("seconds")
	plt.legend(loc='upper left');
	plt.show()
	
	print("h")
	#vec = [np.cos(np.pi/4)*np.sin(np.pi/4),np.sin(np.pi/4)*np.sin(np.pi/4), np.cos(np.pi/4)]
	
	


	#t = spherical_hash_tree(2)2


