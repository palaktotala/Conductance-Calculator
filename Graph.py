import networkx as nx
import matplotlib.pyplot as plt

class Graph:
	def __init__(self,vertices=0):
		self.nodes=vertices
		# self.adjMatrix=[[-1]*self.nodes for i in range(self.nodes)]
		self.adjMatrix=None
		self.shorted=[]
		self.rhs=[]

	def draw_graph_old(self,list_elements):
		g=nx.Graph()
		g.clear()
		for i in list_elements:
			g.add_edge(i[0],i[1],r=i[2])
		pos = nx.spring_layout(g, scale=1.25)
		nx.draw(g,pos,with_labels=True,edge_color='black',width=3,linewidths=3,node_color='red',node_size=1700,font_size=25,font_color="yellow",font_weight="bold")
		edge_labels = nx.get_edge_attributes(g,'r')
		nx.draw_networkx_edge_labels(g, pos, edge_labels = edge_labels,font_size=20,font_color="gray",label_pos=0.3)
		plt.savefig("out_file1.jpg")
		g.clear()
		plt.close()

	def draw_graph_new(self,reduced):
		reduced_elements=reduced
		g=nx.Graph()
		g.clear()
		i=0
		while(i<len(reduced_elements)):
			s=reduced_elements[i][0]
			d=reduced_elements[i][1]
			reduced_elements[i][2]=[reduced_elements[i][2]]
			j=i+1
			while(j<len(reduced_elements)):
				if((reduced_elements[j][0]==s and reduced_elements[j][1]==d) or (reduced_elements[j][1]==s and reduced_elements[j][0]==d)):
					reduced_elements[i][2].append(reduced_elements[j][2])
					del reduced_elements[j]
				j=j+1
			i=i+1
		
		for i in reduced_elements:
			g.add_edge(i[0],i[1],r=i[2])
		pos = nx.spring_layout(g, scale=1.25)
		nx.draw(g,pos,with_labels=True,edge_color='black',width=3,linewidths=3,node_color='red',node_size=1700,font_size=25,font_color="yellow",font_weight="bold")
		edge_labels = nx.get_edge_attributes(g,'r')
		nx.draw_networkx_edge_labels(g, pos, edge_labels = edge_labels,font_size=20,font_color="gray",label_pos=0.3)
		plt.savefig("out_file2.jpg")
		g.clear()
		plt.close()

	def draw_graph_final(self,final_flow):
		print(final_flow)
		i=0
		while i<len(final_flow):
			s=final_flow[i][0]
			d=final_flow[i][1]
			final_flow[i][2]=[abs(round(final_flow[i][2],2))]
			j=i+1
			while(j<len(final_flow)):
				if(final_flow[j][0]==s and final_flow[j][1]==d):
					final_flow[i][2].append(abs(round(final_flow[j][2],2)))
					del final_flow[j]
				j=j+1
			i=i+1
		g=nx.Graph()
		g.clear()
		for i in final_flow:
			g.add_edge(i[0],i[1],c=i[2])
		pos = nx.spring_layout(g, scale=1.25)
		nx.draw(g,pos,with_labels=True,edge_color='blue',width=3,linewidths=3,node_color='red',node_size=1700,font_size=25,font_color="black",font_weight="bold")
		edge_labels = nx.get_edge_attributes(g,'c')
		nx.draw_networkx_edge_labels(g, pos, edge_labels = edge_labels,font_size=20,font_color="gray",label_pos=0.35)
		plt.savefig("out_file3.jpg")
		g.clear()
		plt.close()		

	def setadjMatrix(self):
		self.adjMatrix=[[-1]*self.nodes for i in range(self.nodes)]
		self.rhs=[0.0]*self.nodes

	def add_edge(self,node1,node2,conductance, eps = 1.0/(10**10)):
		if(node1==node2):
			return

		# if (conductance == 0.0):
		# 	conductance = 1.0/eps
		# else:
		# 	conductance = 1.0/float(conductance)

		self.adjMatrix[node1][node2]=conductance
		self.adjMatrix[node2][node1]=conductance

	#changing here		
	def gui_input(self,arrinput):
		for i in range(len(arrinput)):
			if self.adjMatrix[arrinput[i][0]][arrinput[i][1]]==-1:
				self.add_edge(arrinput[i][0],arrinput[i][1],arrinput[i][2])
			else:
				r1=self.adjMatrix[arrinput[i][0]][arrinput[i][1]]
				r2=float(arrinput[i][2])
				equivalent_conductance=(1/r1)+(1/r2)
				# equivalent_conductance= r1 + r2
				
				equivalent_conductance=1/equivalent_conductance
				self.add_edge(arrinput[i][0],arrinput[i][1],equivalent_conductance)

	def make_eqns(self,inlet,outlet,netflow):
		coeffs=[[0.0]*self.nodes for i in range(self.nodes)]

		print(self.adjMatrix)
		for z,i in enumerate(self.adjMatrix):
			
			for j in i:
				if(j!=-1 and coeffs[z][z]!=0.0):
					coeffs[z][z]+=float(1/j)
				elif(j!=-1 and coeffs[z][z]==0.0):
					coeffs[z][z]=float(1/j)

			for index,j in enumerate(i):
				if(index!=z and self.adjMatrix[z][index]!=-1):
					coeffs[z][index]=-float(1/self.adjMatrix[z][index])



		
		self.rhs[inlet]=netflow
		self.rhs[outlet]=-1.0*netflow

		#reducing 1 Variable
		n=max(inlet,outlet)
		#removed variable of Node n by substituting Vn=0.0press
		# newcoeffs=[[0.0]*(self.nodes-1) for i in range(self.nodes)]
		newcoeffs=[]
		for i,c in enumerate(coeffs):
			newcoeffs.append(c[:n]+c[n+1:])

		M = [newcoeffs[i] + [self.rhs[i]] for i in range(len(coeffs))][:-1]
		print(M)
		self.solver(M)

		copy=[x[:] for x in M]
		pressure=self.calc_pressure(copy,max(inlet,outlet))
		

		return abs(M[min(inlet, outlet)][-1]),pressure



	def zero_case(self,list):
		replaced=[]
		for index,i in enumerate(list):
			conductance=i[2]
			if(conductance==0):
				a=min(i[0],i[1])
				b=max(i[0],i[1])
				replaced.append((a,b))
				del list[index]

		self.shorted=sorted(replaced,key=lambda x: x[1], reverse=True)

		no_of_zeroes=len(self.shorted)
		self.nodes=self.nodes-no_of_zeroes#changing total number of nodes

		for (a,b) in self.shorted:
			element=b
			for index,l in enumerate(list):
				if(l[0]==b):
					l[0]=a
				if(l[1]==b):
					l[1]=a
				if(l[0]>b):
					l[0]=l[0]-1
				if(l[1]>b):
					l[1]=l[1]-1

		return list[:]

	def new_inlet(self,inlet):
		
		flag=False
		for (a,b) in self.shorted:
			if(inlet==b):
				inlet=a
				flag=True
				break

		if(not(flag)):
			for index,(a,b) in enumerate(self.shorted):
				if(inlet>b):
					inlet=inlet-len(self.shorted[index:])
					break

		
		return inlet


	def new_outlet(self,out):
		
		flag=False
		for (a,b) in self.shorted:
			if(out==b):
				out=a
				flag=True
				break

		if(not(flag)):
			for index,(a,b) in enumerate(self.shorted):
				if(out>b):
					out=out-len(self.shorted[index:])
					break

		
		return out

	def is_Shorted(self,inlet,outlet):
		return ((inlet,outlet) in self.shorted)

	def solver(self,m, eps = 1.0/(10**10)):
		h=len(m)
		w=len(m[0])
		for y in range(0,h):
			for y2 in range(y+1, h):    # Eliminate column y
				c = m[y2][y] / m[y][y]
				for x in range(y, w):
					m[y2][x] =m[y2][x]- m[y][x] * c
		for y in range(h-1, 0-1, -1): # Backsubstitute
			c  = m[y][y]
			for y2 in range(0,y):
				for x in range(w-1, y-1, -1):
					m[y2][x] = m[y2][x] -  m[y][x] * m[y2][y] / c
			m[y][y] /= c
			for x in range(h, w):       # Normalize row y
				m[y][x] /= c


	def calc_pressure(self,M,x):
		pres=[0]*(len(M)+1)
		j=0
		for i in range(0,len(M)):
			pres[j]=M[i][i]*M[i][len(M)]
			j+=1
			if j==x:
				pres[j]=0
				j+=1
		return pres
      
	def calc_flow_node(self,pres,conductances):
		t=len(pres)
		x=len(conductances)
		flows=[[]for i in range(t)]
		fflo=[[]for i in range(x)]
		power=[[]for i in range(t)]
		fpow=[[]for i in range(x)]
		for i in range(t):
			for j in range(t):
				flows[i].append([])
				power[i].append([])
		for i in conductances:
			pres_diff=pres[i[1]]-pres[i[0]]
			cur=pres_diff/i[2]
			flows[i[0]][i[1]].append(cur)
			powr=cur*cur*i[2]
			power[i[0]][i[1]].append(powr)
		coef=0
		coef1=0
		for i in range(t):
			for j in range(t):
				if flows[i][j]!=[]:
					for k in flows[i][j]:
						fflo[coef].append(i)
						fflo[coef].append(j)
						fflo[coef].append(k)
						
						coef+=1
					for k in power[i][j]:
						fpow[coef1].append(i)
						fpow[coef1].append(j)
						fpow[coef1].append(k)
						coef1+=1


		return fflo,fpow




    









			



			




		














	














	















	
