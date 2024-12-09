import pandas as pd
import numpy as np
from scipy.spatial import KDTree
# The "grid.dat" files having the node info is converted to a "grid.csv" file. #
df = pd.read_csv('grif.csv', header = None) # converting to a dataframe from the "grid.csv" fie
df1 = df    # a copy of the dataframe to add the closest nodes #
nos_elem = df.iloc[0,0]     # number of elements #
nos_nodes = df.iloc[0,1]    # number of nodes #
x = np.zeros(nos_elem)      # x coordinate of each node #
y = np.zeros(nos_elem)      # y coordinate of each node #
z = np.zeros(nos_elem)      # z coordinate of each node #
points = [[[] for _ in range(3)]for _ in range(nos_elem)]   # the coordinates of each node as a 2d array#
for i in range(0,nos_elem):
    x[i] = df.iloc[i+1, 1]
    y[i] = df.iloc[i+1, 2]
    z[i] = df.iloc[i+1, 3]
    points[i][0] = x[i]
    points[i][1] = y[i]
    points[i][2] = z[i]

# print(points)

#       METHOD-1        #
# ------Code for finding the closest nodes (slow and inefficient)----- #

# defining a function to find closest nodes


def closest_nodes(i,n):
     dist = [[0]*nos_elem]*nos_elem
     for j in range(0,nos_elem):
         euclidean = np.sqrt((x[i-1] - x[j])**2 + (y[i-1] - y[j])**2 + (z[i-1] - z[j])**2) # finding euclidean distance
         dist[j] = [j+1, euclidean]
         dist_sorted = sorted(dist, key=lambda x: x[1]) # sorting to find the closest n nodes
         closest_n_nodes = dist_sorted[0:n]
     index_values = [row[0] for row in closest_n_nodes] # returning the index of closest n nodes
     return index_values


for i in range(1, nos_elem):
    for j, value in enumerate(closest_nodes(i,5)):
        df1.iloc[i, j+4] = value
print(df1)

df1.to_csv("gridnew2.dat", sep = ' ', header = None, index = None)
print(closest_nodes(1, 5))
##############################################################################
# The above method is slow and inefficient, using a library in SciPy to find the closest nodes

"""
#       METHOD-2        #
# -----Code for finding the closest nodes (using KDTree)-----#


def closest_nodes(i, k):
    kdtree = KDTree(points)
    distances, indices = kdtree.query(i, k=k)
    indices_list = indices.tolist()     # List of indices of closest "k" nodes to an indice "i" #
    return indices_list


indices_matrix = [[[] for _ in range(5)] for _ in range(nos_elem)]

for i in range(0,nos_elem):
    indices_matrix[i] = [x + 1 for x in closest_nodes(points[i],5)]
    print(i)# a 2D matrix with closest nodes to each node #

# recall: array_2D[row number][column number]

# add to the dataframe

df1 = df

for i in range(0, nos_elem): # range in inclusive
   indices_matrix1 = [int(items) for items in indices_matrix[i]]
   df1.loc[i+1, 4:9] = indices_matrix1 # add to dataframe#

df1 = df1.replace({np.nan: None})
df1 = df1.iloc[0:nos_elem]
df2 = df1
df2.iloc[1:,4:] = df1.iloc[1:,4:].astype(int)
df2.iloc[0,1] = int(df2.iloc[0,1])
df2.to_csv("gridnew1.dat", sep = ' ', header = None, index = None)

"""