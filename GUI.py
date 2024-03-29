import sys
from PyQt5 import QtCore, QtGui, QtWidgets


from Graph import Graph

class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.resize(700, 565)
		MainWindow.setMinimumSize(QtCore.QSize(700, 565))
		MainWindow.setMaximumSize(QtCore.QSize(700, 565))
		self.centralwidget = QtWidgets.QWidget(MainWindow)
		self.centralwidget.setObjectName("centralwidget")
		self.l_nodes = QtWidgets.QLabel(self.centralwidget)
		self.l_nodes.setGeometry(QtCore.QRect(30, 10, 71, 51))
		self.l_nodes.setObjectName("l_nodes")
		self.l_graph = QtWidgets.QLabel(self.centralwidget)
		self.l_graph.setGeometry(QtCore.QRect(30, 70, 121, 41))
		self.l_graph.setObjectName("l_graph")
		self.btn_calculate = QtWidgets.QPushButton(self.centralwidget)
		self.btn_calculate.setGeometry(QtCore.QRect(470, 390, 211, 51))
		self.btn_calculate.setObjectName("btn_calculate")
		self.tv_result = QtWidgets.QLineEdit(self.centralwidget)
		self.tv_result.setGeometry(QtCore.QRect(470, 490, 141, 41))
		self.tv_result.setObjectName("tv_result")
		self.nodes = QtWidgets.QDoubleSpinBox(self.centralwidget)
		self.nodes.setGeometry(QtCore.QRect(100, 20, 81, 31))
		self.nodes.setToolTipDuration(-1)
		self.nodes.setDecimals(0)
		self.nodes.setObjectName("nodes")
		self.line_2 = QtWidgets.QFrame(self.centralwidget)
		self.line_2.setGeometry(QtCore.QRect(10, 450, 441, 20))
		self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
		self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
		self.line_2.setObjectName("line_2")
		self.tv_input = QtWidgets.QLineEdit(self.centralwidget)
		self.tv_input.setGeometry(QtCore.QRect(230, 20, 101, 31))
		self.tv_input.setObjectName("tv_input")
		self.btn_add = QtWidgets.QPushButton(self.centralwidget)
		self.btn_add.setGeometry(QtCore.QRect(350, 20, 91, 31))
		self.btn_add.setObjectName("btn_add")
		self.line_3 = QtWidgets.QFrame(self.centralwidget)
		self.line_3.setGeometry(QtCore.QRect(10, 60, 441, 21))
		self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
		self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
		self.line_3.setObjectName("line_3")
		self.line = QtWidgets.QFrame(self.centralwidget)
		self.line.setGeometry(QtCore.QRect(189, 10, 41, 61))
		self.line.setFrameShape(QtWidgets.QFrame.VLine)
		self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
		self.line.setObjectName("line")
		self.listWidget = QtWidgets.QListWidget(self.centralwidget)
		self.listWidget.setGeometry(QtCore.QRect(460, 40, 231, 331))
		self.listWidget.setViewMode(QtWidgets.QListView.ListMode)
		self.listWidget.setModelColumn(0)
		self.listWidget.setBatchSize(200)
		self.listWidget.setObjectName("listWidget")
		self.line_4 = QtWidgets.QFrame(self.centralwidget)
		self.line_4.setGeometry(QtCore.QRect(410, 10, 81, 581))
		self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
		self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
		self.line_4.setObjectName("line_4")
		self.btn_clear = QtWidgets.QPushButton(self.centralwidget)
		self.btn_clear.setGeometry(QtCore.QRect(310, 390, 121, 51))
		self.btn_clear.setObjectName("btn_clear")
		self.image = QtWidgets.QLabel(self.centralwidget)
		self.image.setEnabled(True)
		self.image.setGeometry(QtCore.QRect(30, 110, 381, 261))
		self.image.setText("")
		self.image.setScaledContents(True)
		self.image.setObjectName("image")
		self.label = QtWidgets.QLabel(self.centralwidget)
		self.label.setGeometry(QtCore.QRect(470, 455, 111, 41))
		self.label.setObjectName("label")
		self.btn_show_old = QtWidgets.QPushButton(self.centralwidget)
		self.btn_show_old.setGeometry(QtCore.QRect(10, 390, 141, 51))
		self.btn_show_old.setObjectName("btn_show_old")
		self.btn_show_new = QtWidgets.QPushButton(self.centralwidget)
		self.btn_show_new.setGeometry(QtCore.QRect(160, 390, 141, 51))
		self.btn_show_new.setObjectName("btn_show_new")
		self.label_2 = QtWidgets.QLabel(self.centralwidget)
		self.label_2.setGeometry(QtCore.QRect(455, 10, 231, 17))
		font = QtGui.QFont()
		font.setPointSize(11)
		self.label_2.setFont(font)
		self.label_2.setObjectName("label_2")
		self.tv_inlet = QtWidgets.QLineEdit(self.centralwidget)
		self.tv_inlet.setGeometry(QtCore.QRect(20, 490, 101, 41))
		self.tv_inlet.setObjectName("tv_inlet")
		self.tv_outlet = QtWidgets.QLineEdit(self.centralwidget)
		self.tv_outlet.setGeometry(QtCore.QRect(150, 490, 91, 41))
		self.tv_outlet.setObjectName("tv_outlet")
		self.label_3 = QtWidgets.QLabel(self.centralwidget)
		self.label_3.setGeometry(QtCore.QRect(30, 470, 67, 17))
		self.label_3.setObjectName("label_3")
		self.label_4 = QtWidgets.QLabel(self.centralwidget)
		self.label_4.setGeometry(QtCore.QRect(150, 470, 81, 17))
		self.label_4.setObjectName("label_4")
		self.line_5 = QtWidgets.QFrame(self.centralwidget)
		self.line_5.setGeometry(QtCore.QRect(450, 370, 261, 20))
		self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
		self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
		self.line_5.setObjectName("line_5")
		self.tv_pressure = QtWidgets.QLineEdit(self.centralwidget)
		self.tv_pressure.setGeometry(QtCore.QRect(270, 490, 121, 41))
		self.tv_pressure.setObjectName("tv_pressure")
		self.label_5 = QtWidgets.QLabel(self.centralwidget)
		self.label_5.setGeometry(QtCore.QRect(270, 470, 121, 21))
		self.label_5.setObjectName("label_5")
		MainWindow.setCentralWidget(self.centralwidget)
		self.statusBar = QtWidgets.QStatusBar(MainWindow)
		self.statusBar.setFocusPolicy(QtCore.Qt.NoFocus)
		self.statusBar.setSizeGripEnabled(True)
		self.statusBar.setObjectName("statusBar")
		MainWindow.setStatusBar(self.statusBar)


		self.list_elements=[]
		self.reduced_elements=[]
		self.graph=Graph()


		self.nodes.valueChanged.connect(self.readNodes)
		self.btn_add.clicked.connect(self.addList)
		self.btn_clear.clicked.connect(self.clearAll)
		self.btn_show_old.clicked.connect(self.drawgraph_old)
		self.btn_show_new.clicked.connect(self.drawgraph_new)
		self.tv_input.returnPressed.connect(self.addList)
		self.btn_calculate.clicked.connect(self.calculate)
		self.listWidget.itemClicked.connect(self.removeitem)



		self.retranslateUi(MainWindow)
		QtCore.QMetaObject.connectSlotsByName(MainWindow)





	def retranslateUi(self, MainWindow):
		_translate = QtCore.QCoreApplication.translate
		MainWindow.setWindowTitle(_translate("MainWindow", "Conductance Calculator"))
		self.l_nodes.setText(_translate("MainWindow", "Nodes"))
		self.l_graph.setText(_translate("MainWindow", "Pipe Network"))
		self.btn_calculate.setToolTip(_translate("MainWindow", "Calculate Net Conductance"))
		self.btn_calculate.setText(_translate("MainWindow", "CALCULATE"))
		self.nodes.setToolTip(_translate("MainWindow", "Enter the Number of Nodes"))
		self.tv_input.setToolTip(_translate("MainWindow", "Inlet, Outlet, Conductance"))
		self.btn_add.setText(_translate("MainWindow", "ADD"))
		self.btn_clear.setToolTip(_translate("MainWindow", "Clear Network and List"))
		self.btn_clear.setText(_translate("MainWindow", "CLEAR"))
		self.image.setToolTip(_translate("MainWindow", "Pipe Network"))
		self.label.setText(_translate("MainWindow", "Net Conductance"))
		self.btn_show_old.setToolTip(_translate("MainWindow", "Display pipe network from given Values"))
		self.btn_show_old.setText(_translate("MainWindow", "OLD NETWORK"))
		self.btn_show_new.setToolTip(_translate("MainWindow", "Display Simplified pipe network"))
		self.btn_show_new.setText(_translate("MainWindow", "NEW NETWORK"))
		self.label_2.setText(_translate("MainWindow", " Inlet, Outlet, Conductance"))
		self.tv_inlet.setToolTip(_translate("MainWindow", "Inlet in the Orignial Network"))
		self.tv_outlet.setToolTip(_translate("MainWindow", "Outlet in The Original Network"))
		self.label_3.setText(_translate("MainWindow", "Inlet"))
		self.label_4.setText(_translate("MainWindow", "Outlet"))
		self.label_5.setText(_translate("MainWindow", "Pressure difference Applied"))
		self.tv_pressure.setToolTip(_translate("MainWindow", "Positive to Vertex"))



	def readNodes(self):
		self.graph.nodes=int(self.nodes.value())
		self.graph.setadjMatrix()


	def removeitem(self):
		index=(self.listWidget.FlowRow())
		self.listWidget.takeItem(index)
		del self.list_elements[index]
		


	def clearAll(self):
		self.listWidget.clear()
		self.nodes.clear()
		self.image.clear()
		self.list_elements=[]
		self.reduced_elements=[]
		self.tv_input.clear()
		self.tv_result.clear()
		self.tv_inlet.clear()
		self.tv_outlet.clear()
		self.tv_pressure.clear()
		self.statusBar.showMessage("")
		self.graph.nodes=0
		self.graph.setadjMatrix()

	def addList(self):
		try:
			self.statusBar.showMessage("")
			str=self.tv_input.text()
			
			self.tv_input.clear()
			# adding to graph
			list=str.split(',')
			list[0]=int(list[0])
			list[1]=int(list[1])
			list[2]=1.0/float(list[2])
			if(list[2]<0):
				raise ValueError
			if(list[0]< 0 or list[0] >=self.graph.nodes or list[1]<0 or list[1]>=self.graph.nodes):
				raise ValueError
			
			self.listWidget.addItem(str)
			self.list_elements.append(list) 
		except(ValueError,IndexError):
			self.statusBar.showMessage("Invalid Input Format")

	def drawgraph_old(self):
  
		copy=[x[:] for x in self.list_elements]
		self.graph.draw_graph_old(copy)

		pixmap=QtGui.QPixmap("out_file1.jpg")
		pixmap=pixmap.scaled(self.image.width(),self.image.height(),QtCore.Qt.KeepAspectRatio)
		self.image.setPixmap(pixmap)

	def drawgraph_new(self):

		copyl=[x[:] for x in self.list_elements]
		self.reduced_elements=self.graph.zero_case(copyl)
		copy=[x[:] for x in self.reduced_elements]
		self.graph.draw_graph_new(copy)
		pixmap=QtGui.QPixmap("out_file2.jpg")
		pixmap=pixmap.scaled(self.image.width(),self.image.height(),QtCore.Qt.KeepAspectRatio)
		self.image.setPixmap(pixmap)

	def drawgraph_final(self,Flows_wire,power_wire):
		
		self.graph.draw_graph_final(Flows_wire)
		pixmap=QtGui.QPixmap("out_file3.jpg")
		pixmap=pixmap.scaled(self.image.width(),self.image.height(),QtCore.Qt.KeepAspectRatio)
		self.image.setPixmap(pixmap)
		for i in power_wire:
			print(i[0]," : ",i[1]," Power Dissipated= ",i[1])


	def calculate(self):
		self.graph.setadjMatrix()
		
		self.graph.gui_input(self.reduced_elements)
		inlet,outlet,pressure=-1,-1,0.0
		
		try:
			inlet=int(self.tv_inlet.text())
			outlet=int(self.tv_outlet.text())
			pressure=float(self.tv_pressure.text())
			if (inlet,outlet) in self.graph.shorted:
				self.tv_result.setText("Zero(Shorted)")
				return
			if(pressure<0.0):
				raise ValueError
			if(inlet < 0 or inlet >=self.graph.nodes):
				raise ValueError
			if(outlet < 0 or outlet >=self.graph.nodes):
				raise ValueError 
		except(ValueError,IndexError,TypeError):
			self.statusBar.showMessage("Invalid Input inlet or outlet")

		ninlet=self.graph.new_inlet(inlet)
		noutlet=self.graph.new_outlet(outlet)
		print("ninlet,noutlet",str(ninlet),str(noutlet))

		net_conductance,pressure_nodes=self.graph.make_eqns(ninlet,noutlet,1)
		netFlow=pressure/net_conductance
		net_pressure,pressure_nodes=self.graph.make_eqns(ninlet,noutlet,netFlow)
		Flows_wire,power_wire=self.graph.calc_flow_node(pressure_nodes,self.reduced_elements)
		print("Flow Supplied=",netFlow)
		# print(Flows_wire)
		# print(power_wire)
		copy=[x[:] for x in Flows_wire]
		self.drawgraph_final(copy,power_wire)
		net_conductance = 1.0/float(net_conductance)
		self.tv_result.setText(str(net_conductance))

def createWindow():
	app=QtWidgets.QApplication(sys.argv)
	MainWindow=QtWidgets.QMainWindow()
	ui=Ui_MainWindow()
	ui.setupUi(MainWindow)
	MainWindow.show()
	sys.exit(app.exec_())

createWindow()
