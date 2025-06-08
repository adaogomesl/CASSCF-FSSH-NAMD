import json
import matplotlib.pyplot as plt 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#def load_data(json_file):
with open('track-dbh-unsub.json','r') as indata:
	data = json.load(indata)

name, cmmt, para, Etot, type, hop, label=data

def states(txt):
	#this function's purpose is to convert a txt file into a  list that will allow you to organize the different trajectories into their respective states
	
	#the following 2 lines create a dictionary, which is used to store the traj information in a more organized way and a counter which
	#allows you to create a key with the proper state
	dic={}
	q = 0
	st_i = list()	
	#the following in the loop opens the txt file that you have create.
	with open(txt, 'r') as infile:
		#the following lines allows you to concacatenate a letter 's' ,representated of 'state' + the number of the state to create the keys for the dictionaries.
		# after it break down txt file by line which allows you to place the traj into their respective keys (states)
		for line in infile:
			#the line below creates the keys
			st = "s"+str(q)
			#the line below splits the txt by line
			dic[st] = line.split()
			#the line below then raises the counter to make the different keys that are representative of their states
			q = q + 1

	#the following loops now breaks down the dictionary into a list which is then returned
	for i in dic: #this is where the loop enters the dictionary through is keys, i = keys, there for to obtaine the list for the specific states is dictionary[key]
		index= dic[i]
		#the following line is an abbreviated form a list, the topic is known as a list comprehension for python, and it allows for the list to be converted from
		#a str value types to int value types
		i = [int(s) for s in index]
		st_i.append(i)
	return st_i

def spaghetti_plt(st_i, cutoff,meas_1,meas_2):
	#this function allows the user to take the list that is returned in the states() function. Here the user is able to input the list, the time length cuttoff, and the measurement
	#that are being plotted. 

	#the following 3 lines are initializing the plots as a local variable rather than a global one. 
	x,y=[],[]
	x2,y2=[],[]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for i in st_i: 
	#this is going into the different states that contain all the traj in that respective state
	#i in this case will be the list of the traj since the inputted list is a 2-D list with the dimesions varying within each state, however the first level.
	#being the number of states. 
	 
		for j in i: 
		#in this loop we are entering each of the lists and i should be the actual number of the trajectory
			j = j - 1 
			#the line above is decreasing one from the trajectory number due to the python starts indexing on 0 rather than 1
			ts = len(para[j])
			
			#the if-else loop below will allow you compare the time length of the trajectory that we are one and determine whether it makes the cutoff
			if ts <= cutoff:
				#if the that specific traj makes the cuttoff then we are able to index into them to get the measurements for every timestep into the x & y lists
				for k in range(ts):
					x.append(para[j][k][0][meas_1][0])
					y.append(para[j][k][0][meas_2][0])
				#once the list has finished going through all the time steps we are able to plot them
				ax.plot(x,y,marker='o',markersize=0.,linewidth=0.5,alpha=0.3)
				x=list()
				y=list()
			else:
				continue


	#this part determines the different hopping points, since we are focus on the hopinh points landing on the S0, thus we are focusing on that state.
	for i in st_i[0]:
		#again we are subtracting one from the actual value within the list since python starts indexing at 0, rather than 1. 
		i = i - 1
		ts = len(para[i])
		#determines whether that specific traj within the s0 list that have hopped and landed on the ground state sastifies the previous cuttoff. 
		if ts <= cutoff:
			#if it does then it index within the hop list to detemine the timesteps which the trajectory hopped in
			ci = hop[i]
			if len(ci)>0: 
			#this determines whether the lists within the hop list (2-dimensional lists) are occupied, if they are then we indexed into them 
				#as per usual we subtract one from the actual timestep that are hopped since the indexing will be off by one. 
				ci = ci[-1]-1
				#with the following lines we are able to determine the actual measurements that are in that timestep, using the same measurements that were in the original
				#spaghetti plot.  
				x2.append(para[i][ci][0][meas_1][0])
				y2.append(para[i][ci][0][meas_2][0])
			else:
				continue

		else:
			continue
	
	#this will actually plot the hopping point on the graph that has the random spaghetti plot as an overlay. It is a seperate from them, and acts as a layer rather than a different graph. 
	ax.scatter(x2,y2,color='black',marker='o',s=3,zorder=3)

	ax.set_xlabel(r'C-N (Å)',color='green',fontsize=14,labelpad=8)
	ax.set_ylabel(r'C-N (Å)',color='blue',fontsize=14,labelpad=8)
	ax.set_xlim(1.3,5.1)
	ax.set_ylim(1.3,5.1)
	plt.savefig("spaghetti-CN-CN.png",dpi=1200)
	###THIS IS WHERE THE X AND Y LIMITS ARE PLACED###
	#ax.set_xlabel(r'Dihedral H-C-C-C (°)',color='blue',fontsize=14,labelpad=8)
	#ax.set_ylabel(r'C-C (Å)',color='red',fontsize=14,labelpad=8)
	#ax.set_xlim(40,210)
	#ax.set_ylim(1.3,2.6)
	#plt.savefig("spaghetti-CC-dihedral.png",dpi=1200)

def avg_sp(st_i, cutoff,avg_1,avg_2,meas_3):
	#this is the same as the spaghetti_plt() function but instead there is the functionality that allows you to average two values within the things 
	#that you are measuring. 
	
	x,y,z,x2,y2,z2,l = [],[],[],[],[],[],[]	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for i in st_i: 
	#this is functioning as the main loop that is indexing within the state lists to pull out the data that we are trying to extract.
		for j in i: 
		#in this loop we are entering each of the lists and j should be the actual number of the trajectory, while i is the list (state) that is being indexed to. 
		#we substract one from the indexed value since the indexing actually starts at zero in python. 
			j = j - 1 
			ts = len(para[j])
			#this is again doing the comparison to determine whether the trajectory's length makes the cutoff or not
			if ts <= cutoff:
				for k in range(ts):
					#this is appends the measurements to the different lists for all the traj's timesteps for all the traj's that made the cutoff.
					x.append(para[j][k][0][avg_1][0])
					y.append(para[j][k][0][avg_2][0])
					z.append(para[j][k][0][meas_3][0])
				
				#the following lines does the zip function, which basically groups the different lists together, which is important for the 
				#list comprehension tool that is used in the following lines. 
				lists = [x,y]
				#the following line converts the tuples that were outputted from the zip function into arrays that allow for array functions
				a = [np.array(x) for x in lists]
				#the following line then uses another list comprehension in order to determine the mean of the two lists for each value within
				#the list, so (value 1 from list 1 + value 1 from list 2)/2
				l = [np.mean(f) for f in zip(*a)]
				print(len(l))	
				#this then plots the average of the two measurements vs the third value that is obtained from the user.
				ax.plot(l,z,marker='o',markersize=0.,linewidth=0.5,alpha=0.3)
				#this then blanks the list so that the list can be used again for the next traj 
				x=list()
				y=list()
				z=list()
				l=list()
				lists=list()
			else:
				continue


	#hopping points here, maybe make this a function 
	for i in st_i[0]:
		i = i - 1
		ts = len(para[i])
		if ts <= cutoff:
			ci = hop[i]
			if len(ci)>0:
				ci = ci[-1] - 1 
				x2.append(para[i][ci][0][avg_1][0])
				y2.append(para[i][ci][0][avg_2][0])
				z2.append(para[i][ci][0][meas_3][0])
			else:
				continue

		else:
			continue

	#this is then doing the same averaging that loops for the entire trajectory, but just for the hopping points.
	#this is necessary to overlay the hopping points onto the plot that has the other traj data.
	lists = [x2,y2]
	a = [np.array(x2) for x2 in lists]
	l = [np.mean(f) for f in zip(*a)]
	ax.scatter(l,z2,color='black',marker='o',s=3,zorder=3)
	
	###THIS IS WHERE THE X AND Y LIMITS ARE PLACED###
	#ax.set_xlim(1.4,2)
	#ax.set_ylim(90,120)


	plt.savefig("avg_sp.png",dpi=1200)

def longest_sp(st_i, cutoff,comp_1,comp_2,meas_3):
	#this is a similar function as before but instead of averaging the measurements, it actually makes a comparison in order to determine
	#which of the bonds is longer. 

	x,y,z,x2,y2,z2,l = [],[],[],[],[],[],[]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for i in st_i: 
	#this is going into into the states, as before in order to index into the 2-D list
		for j in i: 
		#here we are indexing into the individual lists that are actually hold the traj #, so j is the traj # and i is the list of states
		#accounts for the different indexing that is occurring due to the python's indexing starting at zero.
			j = j - 1 
			ts = len(para[j])
			if ts <= cutoff:
			#this if-else statement is accounting to determine whether that specific traj that we are one is actually hitting the cutoff. 
				for k in range(ts):
				#if it is then this for-loop iterates through the length of the trajectory and appends the measurements at each timestep to the
				#different lists. 
					x.append(para[j][k][0][comp_1][0])
					y.append(para[j][k][0][comp_2][0])
					z.append(para[j][k][0][meas_3][0])
				for i in range(len(x)):
				#this is the comparison part that is determining which of the bonds is longer, to do this we compare the x and y lists and then append
				#a new list with the longest of the both. This method ensures the correct dimensionalities. 
					if x[i] > y[i]:
						l.append(x[i])
					else:
						l.append(y[i])

				#this then plots the data that has been determined with the other functions, the longest carbons in the l list and then the third measurement. 
				ax.plot(l,z,marker='o',markersize=0.,linewidth=0.5,alpha=0.3)

				# the next lines then clear the list for the next iteration of the code, so for the next traj. 
				x=list()
				y=list()
				z=list()
				l=list()
			else:
				continue


	#hopping points here
	for i in st_i[0]:
		i = i - 1
		ts = len(para[i])
		if ts <= cutoff:
			ci = hop[i]
			if len(ci)>0:
				ci = ci[-1] - 1 
				x2.append(para[i][ci][0][comp_1][0])
				y2.append(para[i][ci][0][comp_2][0])
				z2.append(para[i][ci][0][meas_3][0])
			else:
				continue

		else:
			continue

	#determines which was the longest carbon bond for the hopping point, which takes two lists x2 and y2, and compares them to append a new one l
	for i in range(len(x2)):
		if x2[i] > y2[i]:
			l.append(x2[i])
		else: 
			l.append(y2[i])

	#plots the surface hopping points on top of the previous plots that contain all the other trajectory data. 
	ax.scatter(l,z2,color='black',marker='o',s=3,zorder=3)

	###THIS IS WHERE THE X AND Y LIMITS ARE PLACED###
	#longest C-N versus C-C
	#ax.set_xlim(1.3,22)
	#ax.set_ylim(1.25,2.75)
	#ax.set_xlabel(r'C-N (Å)',color='blue',fontsize=14,labelpad=8)
	#ax.set_ylabel(r'C-C (Å)',color='red',fontsize=14,labelpad=8)

	#longest C-N versus angle C-C-C
	#ax.set_xlim(1.2,25)
	#ax.set_ylim(50,120)
	#ax.set_xlabel(r'Dihedral H-C-C-C (°)',color='blue',fontsize=14,labelpad=8)
	#ax.set_ylabel(r'C-C (Å)',color='red',fontsize=14,labelpad=8)
	plt.savefig("CNxCN.png",dpi=1200)

print(cmmt)
st_i = states('states.txt')
spaghetti_plt(st_i,2000,1,2)
#spaghetti_plt(st_i,2000,5,0)

