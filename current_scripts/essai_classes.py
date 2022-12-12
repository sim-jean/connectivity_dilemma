# QItem for current location and distance
import ast
import statistics as stat
import itertools

import pandas as pd

import current_scripts.functions_project as functions_project
exec(open('current_scripts/functions_project.py').read())

import networkx as nx
import numpy as np

import current_scripts.functions_project as functions_project
import seaborn as sns
import matplotlib.pylab as plt
# from source location
class QItem:
    def __init__(self, row, col, dist):
        self.row = row
        self.col = col
        self.dist = dist

    def __repr__(self):
        return f"QItem({self.row}, {self.col}, {self.dist})"


def minDistance(grid,sourcex,sourcey,destx,desty):
    source = QItem(sourcex, sourcey, 0)


    # To maintain location visit status
    visited = [[False for _ in range(len(grid[0]))]
               for _ in range(len(grid))]

    # applying BFS on matrix cells starting from source
    queue = []
    queue.append(source)
    visited[source.row][source.col] = True
    while len(queue) != 0:
        source = queue.pop(0)

        # Destination found;
        if(source.row==destx and source.col==desty):
            return( source.dist)
            #print("finito")

        # moving up left
        if isValid(source.row - 1, source.col-1, grid, visited):
            queue.append(QItem(source.row - 1, source.col-1, source.dist + 1))
            visited[source.row - 1][source.col-1] = True

        # moving up right
        if isValid(source.row - 1, source.col+1, grid, visited):
            queue.append(QItem(source.row - 1, source.col+1, source.dist + 1))
            visited[source.row - 1][source.col+1] = True

        # moving down right
        if isValid(source.row + 1, source.col + 1, grid, visited):
            queue.append(QItem(source.row + 1, source.col + 1, source.dist + 1))
            visited[source.row + 1][source.col + 1] = True
        # moving down left
        if isValid(source.row + 1, source.col -1, grid, visited):
            queue.append(QItem(source.row - 1, source.col + 1, source.dist + 1))
            visited[source.row + 1][source.col - 1] = True

        # moving up
        if isValid(source.row - 1, source.col, grid, visited):
            queue.append(QItem(source.row - 1, source.col, source.dist + 1))
            visited[source.row - 1][source.col] = True

        # moving down
        if isValid(source.row + 1, source.col, grid, visited):
            queue.append(QItem(source.row + 1, source.col, source.dist + 1))
            visited[source.row + 1][source.col] = True

        # moving left
        if isValid(source.row, source.col - 1, grid, visited):
            queue.append(QItem(source.row, source.col - 1, source.dist + 1))
            visited[source.row][source.col - 1] = True

        # moving right
        if isValid(source.row, source.col + 1, grid, visited):
            queue.append(QItem(source.row, source.col + 1, source.dist + 1))
            visited[source.row][source.col + 1] = True

    return 1000
    #dealt with 0 if no path exists, but it is problematic ; if all patches have the same size
    # and no

# checking where move is valid or not
def isValid(x, y, grid, visited):
    if ((x >= 0 and y >= 0) and
            (x < len(grid) and y < len(grid[0])) and
            (grid[x][y] != 0) and (visited[x][y] == False)):
        return True
    return False





class Land:
    def __init__(self,shape):
        self.shape= shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = functions_project.to_matrix(shape,size)
        self.biod = (shape>=1).dot(functions_project.connectivity_mat(size)).dot(np.transpose((shape>=1)))-sum(shape>=1)
        self.fuel = (shape>=2).dot(functions_project.connectivity_mat(size)).dot(np.transpose((shape>=2)))-sum(shape>=2)
        self.nonz = sum(shape>0)
        self.nb2 = sum(shape==2)
        self.zero = sum(shape==0)
        self.habitat_area_ratio = sum(self.shape>=1)/(self.size**2)


    def view(self):
        ax = plt.axes()
        ax.set_title('Initial biodiv = '+str(self.biod))
        sns.heatmap(self.mat,ax = ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
        plt.show()
    def view_succession(self,index):
        if len(index)==1:
            ax = plt.axes()
            ax.set_title('Biodiv constraint = ' + str(index[0])+" with budget = "+str(budget))
            y = functions_project.to_matrix(self.succession_land[index[0]], self.size)
            sns.heatmap(y, ax=ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
            plt.show()
        else:
            number_x = 0
            rows = len(index)//2
            for i in range(rows):
                for j in range(2):
                    ax = plt.subplot2grid((rows, 2), (i, j))
                    ax.title.set_text("Biodiv constraint =" + str(index[number_x]) + ' with budget ='+str(budget))
                    y = functions_project.to_matrix(self.succession_land[index[number_x]], self.size)
                    ax = sns.heatmap(y, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
                    number_x += 1
            plt.show()
    def equivalent(self):
        m2 = np.rot90(self.mat)
        m3 = np.rot90(m2)
        m4 = np.rot90(m3)
        if self.size==5:
             # Horizontal symmetry
            hor_sym = [self.mat[4], self.mat[3], self.mat[2], self.mat[1], self.mat[0]]
        elif self.size==4:
            hor_sym=[self.mat[3], self.mat[2], self.mat[1], self.mat[0]]
        elif self.size==3:
            hor_sym = [self.mat[2], self.mat[1], self.mat[0]]
            # 3 clockwise transformations
        else:
            ValueError("Size is not supported")
        m6 = np.rot90(hor_sym)
        m7 = np.rot90(m6)
        m8 = np.rot90(m7)
        evalu = (tuple(list(itertools.chain(*m2))),
                tuple(list(itertools.chain(*m3))),
                tuple(list(itertools.chain(*m4))),
                tuple(list(itertools.chain(*hor_sym))),
                tuple(list(itertools.chain(*m6))),
                tuple(list(itertools.chain(*m7))),
                tuple(list(itertools.chain(*m8))))
        return set(evalu)

    def fuel_dynamics(self,presc_burn):
        post = (np.array(self.shape) + 1) * (1 - np.array(presc_burn))
        land_post = np.minimum(post, 2)
        self.fuel_dynamics=land_post

    def succession(self,constraint_biodiv):
        if len(self.fuel_dynamics)>1:
            biod =np.mat(np.diagonal((self.fuel_dynamics >= 1).dot(functions_project.connectivity_mat(self.size)).dot(np.transpose((self.fuel_dynamics >= 1))) - (self.fuel_dynamics >= 1).sum(
                    axis=1))).T
            fuel = np.transpose(np.matrix(np.diagonal(
                (self.fuel_dynamics >= 2).dot(connectivity_mat(self.size)).dot(np.transpose((self.fuel_dynamics >= 2))) - (self.fuel_dynamics >= 2).sum(axis=1))))

            value = fuel + self.fuel + 1000*(biod<=constraint_biodiv)
            a = value.min(0)
            b = value.argmin(0)

            self.succession_land = {k:v for (k,v) in zip(constraint_biodiv, [np.array(x) for x in self.fuel_dynamics[np.asarray(b, int)][0].tolist()])}
            self.succession_value ={k:v for (k,v) in zip(constraint_biodiv, np.asarray(a)[0].tolist())}
            return self.succession_land, self.succession_value

        else :
            print("prout")

    def shortest_path(self,sourcex,sourcey,destx,desty):
        grid = (self.mat>=1)
        grid = grid.astype(int)
        a = minDistance(grid,sourcex,sourcey,destx,desty)
        b = minDistance(grid,destx,desty,sourcey,sourcex)
        if a!=b:
            return(min(a,b))
        else:
            return(a)


    def IIC_fuel(self):
        paths = self.test_algo(2)
        invert = [1/(1+x) for x in list(paths.values())]
        iic = sum(invert)/(self.size)**4
        return(iic)
    def IIC_biod(self):
        paths = self.test_algo(1)

        invert = [1/(1+x) for x in list(paths.values())]
        iic = sum(invert)/(self.size)**4
        return(iic)
    def characteristic_path_length_biod(self):
        paths = self.test_algo(1)
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        return stats.mean(paths2)
    def characteristic_path_length_fire(self):
        paths = self.test_algo(2)
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]


        return stats.mean(paths2)


    def test_algo(self,k):

        grid = (self.mat >= k)
        grid = grid.astype(int)

        all_coordinates = []
        for i in range(4):
            all_coordinates.append(list(zip([i] * 4, range(4))))
        all_coordinates = [item for sublist in all_coordinates for item in sublist]

        paths = {}
        for i in all_coordinates:
            for j in all_coordinates:
                paths[(i,j)]=minDistance(grid, i[0], i[1], j[0], j[1])

        b = [key for key in list(paths.keys()) if paths[key] == 1000]

        for key in b:
            if paths[(key[1], key[0])] != 1000:
                paths[(key[0], key[1])] = paths[(key[1], key[0])]
        return(paths)

    def diameter_biod(self):
        paths = self.test_algo(1)

        paths2 = [paths[x] for x in paths.keys() if paths[x]<1000]
        return max(paths2)
    def diameter_fire(self):
        paths = self.test_algo(2)

        paths2 = [paths[x] for x in paths.keys() if paths[x]<1000]
        return max(paths2)

    def diagonal_biod(self,direction):
        paths = self.test_algo(1)
        if direction == 'NO-SE':
            return(paths[((0,0),(self.size-1,self.size-1))])
        elif direction=='SO-NE':
            return(paths[((size-1,0),(0,size-1))])

    def adjacency_biod(self):
        paths = self.test_algo(1)


costs = functions_project.cost(size)
budget=1
size=4

pot_fire_value  = (0,1)
pot_fire        = list(itertools.product(pot_fire_value,repeat=size**2))
pot_fire_budget = [x for x in pot_fire if sum(sum(to_matrix(x,size)*costs))<=budget]

presc_burn = pot_fire_budget

landscape = Land(shape=np.array((1,0,2,0,1,2,0,1,2,2,0,2,0,1,2,0)))
landscape.fuel_dynamics(presc_burn)

landscape.succession([2,4,12,24,48])
landscape.view()
landscape.view_succession([2])

# Voir si on peut faire la même chose qu'avec les autres algo en opti notamment
landscape.IIC_biod()
landscape2 = Land(shape=np.array((1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)))
landscape2.IIC_biod()
landscape3 = Land(shape=np.array((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
landscape3.IIC_biod()
landscape4 = Land(shape=np.array((1,1,1,1,1,0,0,1,1,0,0,1,1,1,1,1)))

#essai avec networkx
#find perimeter of habitat in

landscape9 = Land(shape = np.array((1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1)))

def number_empty_neighbors(mat,i,j,size=3,var="biod"):
    if var=="biod":
        threshold=1
    elif var=="fire":
        threshold=2
    if i>0:
        up = (mat[i-1][j]<threshold)
    elif i==0:
        up=1

    if i<size:
        down = (mat[i+1][j]<threshold)
    elif i==size:
        down=1
    if j>0:
        left = (mat[i][j-1]<threshold)
    elif j==0:
        left=1
    if j<size:
        right= (mat[i][j+1]<threshold)
    elif j==size:
        right=1
    return(sum([up,down,right,left]))

def perimeter_biod(mat):
    perimeter = 0
    for j in range(len(mat[0])):
        for i in range(len(mat[0])):
            perimeter+= (mat[i][j]>0)*number_empty_neighbors(mat,i,j)
    return perimeter

def perimeter_fire(mat):
    perimeter = 0
    for j in range(len(mat[0])):
        for i in range(len(mat[0])):
            perimeter+= (mat[i][j]>1)*number_empty_neighbors(mat,i,j, var='fire')
    return perimeter

def landscape_shape_index(mat):
    return perimeter(mat)/16

def perimeter_to_area(mat):
    return perimeter(mat)/sum(sum(mat==1))

#
data2 = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/successions_matched/land4_budget_1_10_8_cut_0_biodiv_index_2.csv")

successions = data2.iloc[0,1:21]
successions = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]

new_dat = {}
new_dat['land']= [successions[x].shape for x in range(len(successions))]
new_dat['perim_biodiv']=[perimeter_biod(successions[x].mat) for x in range(len(successions))]
new_dat['biod_indic']= [successions[x].biod for x in range(len(successions))]
new_dat['fuel_indic']= [successions[x].fuel for x in range(len(successions))]
new_dat['perim_fire']=[perimeter_fire(successions[x].mat) for x in range(len(successions))]
new_dat['charact_path_biod']=[successions[x].characteristic_path_length_biod() for x in range(len(successions))]
new_dat['charact_path_fire']=[successions[x].characteristic_path_length_fire() for x in range(len(successions))]
new_dat['diameter_biod']=[successions[x].diameter_biod() for x in range(len(successions))]
new_dat['diameter_fire']=[successions[x].diameter_fire() for x in range(len(successions))]
new_dat['IIC_biod']=[successions[x].IIC_biod() for x in range(len(successions))]
new_dat['IIC_fuel']=[successions[x].IIC_fuel() for x in range(len(successions))]


data_eval = pd.DataFrame(new_dat)
landscape9.II