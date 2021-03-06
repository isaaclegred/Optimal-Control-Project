import numpy as np
from Utilities import *
from matplotlib import pyplot as plt
import networkx
x = np.linspace(-1,1,10)
y = np.linspace(-1,1,10)
X,Y = np.meshgrid(x, y)
positions = np.stack([X.ravel(), Y.ravel()])
target_graph = Graph_from_Adjacency(np.diag(np.ones(99), 1) + np.diag(np.ones(99), -1)
                                    + np.diag(np.ones(90), 10) + np.diag(np.ones(90),
                                                                        -10) + \
                                    np.diag(np.ones(98), 2) + np.diag(np.ones(98), -2)
                                    + np.diag(np.ones(95), 5) + np.diag(np.ones(95),
                                                                         -5),
                                    positions)

G= networkx.from_numpy_matrix(np.diag(np.ones(99), 1) + np.diag(np.ones(99), -1)
                                    + np.diag(np.ones(90), 10) + np.diag(np.ones(90),
                                                                        -10))
nxpos = {}
for integer in range(100):
    nxpos[integer] = positions[:, integer]

networkx.draw(G, nxpos)

plt.draw()
plt.figure()
target = target_graph[0]
print("target is", target)
for node  in target.N:
    print("neighbor has neighbors", node.N)
LC = LeastCost(lambda x: (2-np.linalg.norm(x)**2), target)
positions, costs = LC.Find_Values(target_graph)
X = positions[:,0]
Y = positions[:, 1]
print(len(costs))
print(len(X))
print("target position is", target.X)

plt.scatter(X,Y, c = costs),
plt.show()
#plt.contourf(np.reshape(costs, (10,10)) (X, Y))
