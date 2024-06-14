import numpy as np
import pylab as pl

class Truss:
    
    def __init__(self):
        self.nodes = {}
        self.connections = {}
        self.loads = {}
        
    # sets a node, adds it to the dictionary if not already exists
    def setNode(self, node_id, x, y, isFixed):
        self.nodes[node_id] = (x,y,isFixed)
        self.loads[node_id] = (0.0,0.0)
        #print(self.nodes)
    
    def setLoad(self, node_id, xf, yf):
        self.loads[node_id] = (xf,yf)
        
    # removes a node and any associated connections
    def removeNode(self, node_id):
        del self.nodes[node_id]
        del self.loads[node_id]
        
        l = list(self.connections.keys())
        for key in l:
            if node_id in key:
                del self.connections[key]
                

        
    # adds a new connection to the current dictionary
    def addConnection(self, firstNode, secondNode):
        if firstNode > secondNode:
            self.connections[(secondNode, firstNode)] = 0

            return
        self.connections[(firstNode,secondNode)] = 0

    
    # removes a connection to the current dictionary
    def removeConnection(self, firstNode, secondNode):
        if firstNode > secondNode:
            del self.connections[(secondNode, firstNode)]
            return
        del self.connections[(firstNode,secondNode)]
        
    def getForces(self):
        self.numForces = len(self.connections)
        self.nonFixedNodes = {}
        
        for key in list(self.nodes.keys()):
            if (not self.nodes[key][2]):
                self.nonFixedNodes[key] = self.nodes[key]
                
        self.numUnknownNodes = len(self.nonFixedNodes)
        
        self.forceIndices = list(self.connections.keys())
        self.nodeIndices = list(self.nonFixedNodes.keys())
        
        # set up the way the force directions affect each point
        self.directions = np.zeros((self.numForces,len(self.nonFixedNodes)),dtype=float).T
        
        for i in range(len(self.nonFixedNodes)):
            curr_id = self.nodeIndices[i]
            
            for j in range(self.numForces):
                if(curr_id in self.forceIndices[j]):
                    if(curr_id == self.forceIndices[j][0]):
                        self.directions[i][j] = -1
                    else:
                        self.directions[i][j] = 1
                else:
                    self.directions[i][j] = 0
        self.lengths = np.zeros(self.numForces, dtype=float)
        
        self.components = np.zeros((self.numForces, 2), dtype=float)
        for i in range(self.numForces):
            member = self.forceIndices[i]

            self.lengths[i] = np.sqrt((self.nodes[member[0]][0] - self.nodes[member[1]][0])**2.0 
                                      + (self.nodes[member[0]][1] - self.nodes[member[1]][1])**2.0)
            
            self.components[i][0] = (self.nodes[member[0]][0] - self.nodes[member[1]][0]) / self.lengths[i]
            self.components[i][1] = (self.nodes[member[0]][1] - self.nodes[member[1]][1]) / self.lengths[i]
        
        
        # construct [A] matrix
        
        A = np.zeros((self.numForces,2*self.numUnknownNodes),dtype=float).T
        
        for i in range(2*self.numUnknownNodes):
            for j in range(self.numForces):
                A[i][j] = self.components[j][int(i) % 2] * self.directions[int(i / 2)][j]
                
        A *= -1.0
        
        k = np.zeros((self.numForces,self.numForces), dtype=float)
        
        for i in range(self.numForces):
            k[i][i] = 1.0 / self.lengths[i]
            
        K = A @ k @ A.T
        
        X = np.zeros(2*self.numUnknownNodes,dtype=float).T
        
        for i in range(self.numUnknownNodes):
            X[2*i] = self.loads[self.nodeIndices[i]][0]
            X[2*i + 1] = self.loads[self.nodeIndices[i]][1]
        U = np.linalg.solve(K,X)
        
        Delta = A.T @ U
        F = k @ Delta
        
        return F

mytruss = Truss()

mytruss.setNode(0,0,0,True)
mytruss.setNode(1,1,2,False)
mytruss.setNode(2,2,0,True)
mytruss.setNode(3,3,0,True)
mytruss.setLoad(1,0,-100)
mytruss.addConnection(0,1); mytruss.addConnection(1,2);
mytruss.addConnection(3, 1)
print(mytruss.getForces())