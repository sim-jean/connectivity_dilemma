import ast
import os

import current_scripts.functions_project as functions_project
exec(open('current_scripts/functions_project.py').read())

class Land_lite:
    def __init__(self,shape):
        self.shape= shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = to_matrix(shape,size)
        self.biod = (shape>=1).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=1)))-sum(shape>=1)
        #self.biod_relative = self.biod/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape))))-sum(np.repeat(1,len(self.shape))))
        self.fuel = (shape>=2).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=2)))-sum(shape>=2)
        #self.fuel_relative = self.fuel/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape))))-sum(np.repeat(1,len(self.shape))))
        self.node_biod = sum(shape>0)
        self.node_fuel = sum(shape==2)
        self.zero = sum(shape==0)
        self.habitat_area_ratio = sum(self.shape>=1)/(self.size**2)
        self.fuel_area_ratio = sum(self.shape==2)/(self.size**2)
        self.test_algo_0 = self.test_algo(0)
        self.test_algo_1 = self.test_algo(1)
        self.test_algo_2 = self.test_algo(2)
        self.components_biod_tot = self.components()
        self.components_fuel_tot = self.components('fuel')
    # Functions for characterization
    def test_algo(self,k):
        if k>=0:
            grid = (self.mat >= k)
        elif k==0:
            grid = (self.mat==k)
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
    def shortest_path(self,sourcex,sourcey,destx,desty):
        grid = (self.mat>=1)
        grid = grid.astype(int)
        a = minDistance(grid,sourcex,sourcey,destx,desty)
        b = minDistance(grid,destx,desty,sourcey,sourcex)
        if a!=b:
            return(min(a,b))
        else:
            return(a)

    def IIC_with_nc(self,var="biod"):
        if var=="biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        invert = [1 / (1 + x) for x in list(paths.values())]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size) ** 4
        return(iic)
    def IIC_without_nc(self,var="biod"):
        if var=="biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        screen1 = [x for x in list(paths.values()) if x<1000]
        invert = [1 / (1 + x) for x in screen1]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size) ** 4
        return(iic)
    def CPL(self,var="biod"):
        if var == "biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        return stats.mean(paths2)
    def CPL_without_nc(self,var="biod"):
        if var == "biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        paths2 = [x for x in paths2 if x<1000]
        try:
            result = stats.mean(paths2)
        except statistics.StatisticsError:
            result = float("nan")
        return result

    def diameter(self,var="biod"):
        if var=="zero":
            paths = self.test_algo_0
        if var=="biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        paths2 = [paths[x] for x in paths.keys() if paths[x]<1000]
        return max(paths2)


    def landscape_shape_index(self,var="biod"):
        if var=='biod':
            return 0.25*self.perimeter(var="biod")/self.size
        elif var=='fuel':
            return 0.25*self.perimeter(var="fuel")/self.size
    def landscape_perimeter_area(self,var="biod"):
        if var=='biod':
            return 0.25*self.perimeter(var="biod")/sum(self.shape>=1)
        elif var=='fuel':
            return 0.25*self.perimeter(var="fuel")/sum(self.shape==2)

    # Description
    # Utility functions : empty neighbors and zero adjacent
    def number_empty_neighbors(self,i,j,var="biod"):
        """
        4 direction empty neighbors
        :param i: Coordinate in grid
        :param j: Coordinate in grid
        :param var: "biod" or "fuel" depending on which metric is evaluated
        :return: number of empty 4 neighbors
        """
        size = self.size - 1
        if var=="biod":
            threshold=1
        elif var=="fuel":
            threshold=2

        if i>0:
            up = (self.mat[i-1][j]<threshold)
        elif i==0:
            up=1

        if i<size:
            down = (self.mat[i+1][j]<threshold)
        elif i==size:
            down=1

        if j>0:
            left = (self.mat[i][j-1]<threshold)
        elif j==0:
            left=1

        if j<size:
            right= (self.mat[i][j+1]<threshold)
        elif j==size:
            right=1

        return(sum([up,down,right,left]))
    def number_zero_adjacent(self,i,j, var="biod"):
        """
        Number of 8 neighbors empty neighbors of cell (i,j) in grid
        :param i: Coordinate in grid
        :param j: Coordinate in grid
        :param var: "biod" or "fuel" depending on which metric is evaluated
        :return:
        """
        size = self.size - 1

        if var=="biod":
            threshold=1
        elif var=="fuel":
            threshold=2

        if i<0 or j<0 or i>size or j>size:
            return(0)

        if i>0:
            up = int(self.mat[i-1][j]==0)
        elif i<=0:
            up=1

        if i<size:
            down = int(self.mat[i+1][j]==0)
        elif i>=size:
            down=1
        if j>0:
            left = int(self.mat[i][j-1]==0)
        elif j<=0:
            left=1
        if j<size:
            right= int(self.mat[i][j+1]==0)
        elif j>=size:
            right=1

        if i>0 and j<size:
            up_right = int(self.mat[i-1][j+1]==0)
        else:
            up_right=1

        if i>0 and j>0:
            up_left = int(self.mat[i-1][j-1]==0)
        else:
            up_left =1

        if i<size and j>0:
            down_left = int(self.mat[i+1][j-1]==0)
        else:
            down_left=1

        if i<size and j<size:
            down_right = int(self.mat[i+1][j+1]==0)
        else:
            down_right=1

        sum_neighbor = sum([up,down,right,left,up_right,up_left,down_right,down_left])
        return(sum_neighbor*(self.mat[i][j]>=threshold))
    def nb_patches(self,var="biod"):
        size=self.size
        nb = []
        if var=="biod":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i,j,var="biod"))
        elif var=="fuel":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i,j,var="fuel"))

        return(1+sum([x==8 for x in nb]))
        # Components
    def components(self,var = "biod"):
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        return g.countIslands()

    def components_perimeter(self, var= "biod",components_var=False):
        vars = var
        if var == "biod":
            first_nodes = self.components_biod_tot[2]
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            first_nodes = self.components_fuel_tot[2]
            check = (self.mat==2).astype(int)
        if len(first_nodes) <= 1:
            if var=="biod":
                return(self.perimeter(var="biod"))
            elif var=='fuel':
                return(self.perimeter(var="fuel"))
        else:
            components = []
            perimeter_components = []
            for nodes in first_nodes:
                queue = []
                queue.append(nodes)
                visited = np.ones((4, 4), dtype=bool) == False
                patches = []
                cols = [0, 0, -1, 1, -1, -1, 1, 1]
                rows = [-1, 1, 0, 0, 1, -1, 1, -1]
                while len(queue) > 0:
                    node = queue.pop(0)

                    for i, j in list(zip(rows, cols)):
                        if isValid(node[0] + i, node[1] + j, check, visited):
                            queue.append((node[0] + i, node[1] + j))
                            visited[node[0] + i, node[1] + j] = True
                            patches.append((node[0] + i, node[1] + j))
                components.append(patches)
                perimeter = 0
                for patch in patches:
                    perimeter += check[patch[0], patch[1]]* self.number_empty_neighbors(patch[0], patch[1], var=vars)

                perimeter_components.append(perimeter)
            if components_var== True:
                return(perimeter_components,components)
            else:
                return (perimeter_components)
    def components_shape_index(self,var="biod"):
        if var =="biod":
            if self.components_biod_tot[0]>1:
                area = max(self.components_biod_tot[1])
                candidate = self.components_biod_tot[1].index(area)
                try:
                    perimeter = self.components_perimeter(var="biod")[candidate]
                except IndexError:
                    perimeter =  self.components_perimeter(var="biod")
            else:
                area = self.components_biod_tot[1][0]
                perimeter = self.components_perimeter(var="biod")

        elif var=="fuel":
            if self.components_fuel_tot[0]>1:
                area = max(self.components_fuel_tot[1])
                candidate = self.components_fuel_tot[1].index(area)
                try:
                    perimeter = self.components_perimeter(var="fuel")[candidate]
                except IndexError:
                    perimeter =  self.components_perimeter(var="fuel")
            else:
                area = self.components_fuel_tot[1][0]
                perimeter = self.components_perimeter(var="fuel")
        return(0.25*perimeter/math.sqrt(area))

        # Overall graph

    def node_degree(self,i,j,var="biod"):
        if i<0 or j<0 or j>self.size-1 or i>self.size-1:
            return(0)
        if var=="biod":
            if self.mat[i,j]==0:
                return(0)
            else:
                return 8- self.number_zero_adjacent(i,j,var='biod')
        elif var=="fuel":
            if self.mat[i,j]==0:
                return(0)
            else:
                return 8- self.number_zero_adjacent(i,j,var='fuel')
    def connectivity_correlation(self,var="biod",option_plot="No"):
        store_node = []
        store_neighbors = []
        for i in range(self.size):
            for j in range(self.size):
                first_coord = [-1,1,0,0,-1,-1,1,1]
                second_coord = [0,0,-1,1,-1,1,-1,1]
                store = []
                if var=="biod":
                    for nb in range(8):
                        store.append(self.node_degree(i+first_coord[nb],j+second_coord[nb],var="biod"))
                    store_node.append(self.node_degree(i,j,var="biod"))
                elif var=="fuel":
                    for nb in range(8):
                        store.append(self.node_degree(i+first_coord[nb],j+second_coord[nb],var="fuel"))
                    store_node.append(self.node_degree(i,j,var="fuel"))

                positive_x=[x for x in store if x!=0]
                try:
                    store_neighbors.append(sum(positive_x)/len(positive_x))
                except stats.StatisticsError:
                    store_neighbors.append(0)
                except ZeroDivisionError:
                    store_neighbors.append(0)

        try:
            coefficient = np.corrcoef(np.array(store_node),np.array(store_neighbors))
        except RuntimeWarning:
            coefficient = float("nan")
        if option_plot=="Yes":
            plt.scatter(store_node,store_neighbors)
            plt.ylabel("Average # of edges for 8 neighbors ")
            plt.xlabel("Number of edges of node")
            plt.show()
        else:
            return coefficient[0,1]
    def perimeter(self,var="biod"):
        perimeter = 0
        if var=="biod":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter+= (self.mat[i][j]>0)*self.number_empty_neighbors(i,j)
        elif var=="fuel":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter += (self.mat[i][j] > 1) * self.number_empty_neighbors(i, j, var='fuel')
        return perimeter

    def perimeter_to_area(self,var='biod'):
        if var=='biod':
            return self.perimeter(var="biod")/math.sqrt(sum(sum(self.mat==1)))
        elif var=='fuel':
            return self.perimeter(var="fuel")/math.sqrt(sum(sum(self.mat==2)))

times = []
timess=0
while timess<10:
    data = pd.read_csv(
        "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_2/successions_matched/land4_budget_2_9_5_cut_9_biodiv_index_12.csv").drop(
        columns=['Unnamed: 0'])

    data = data.iloc[:, 0:20]
    names = list(range(20))
    data.columns = names

    data["index"] = list(range(len(data)))
    # data['index_biod'] = [12]*len(data)
    # could be better to do in the end, with initial numbers as well.

    checker = data.melt(id_vars=["index"])
    checker["variable"] = checker['variable'].astype(int)
    checker = checker.sort_values(["index", "variable"])
    checker = checker.reset_index(drop=True)

    values = list(checker.value)
    values_unique = list(set(values))

    # region Setting up of data storage dictionnaries
    nodes_biod = {}
    nodes_fuel = {}

    biod_score = {}
    fuel_score = {}

    land_perimeter_biod = {}
    land_perimeter_fuel = {}

    land_habitat_area_ratio = {}
    land_fuel_area_ratio = {}

    land_shape_index_biod = {}
    land_shape_index_fuel = {}

    land_diameter_biod = {}
    land_diameter_fuel = {}

    land_connectivity_correlation_biod = {}
    land_connectivity_correlation_fuel = {}

    land_IIC_nc_biod = {}
    land_IIC_nc_fuel = {}

    land_IIC_biod = {}
    land_IIC_fuel = {}

    land_CPL_nc_biod = {}
    land_CPL_nc_fuel = {}

    land_CPL_biod = {}
    land_CPL_fuel = {}

    components_number_biod = {}
    components_number_fuel = {}

    components_area_biod_max = {}
    components_area_fuel_max = {}

    components_shape_index_biod = {}
    components_shape_index_fuel = {}
    # endregion

    # region Compute statistics
    debut = time.time()
    for shape in values_unique:
        land = Land_lite(np.array(ast.literal_eval(shape)))

        nodes_biod[shape] = land.node_biod
        nodes_fuel[shape] = land.node_fuel

        biod_score[shape] = land.biod
        fuel_score[shape] = land.fuel

        land_perimeter_biod[shape] = land.perimeter()
        land_perimeter_fuel[shape] = land.perimeter(var="fuel")

        land_habitat_area_ratio[shape] = land.habitat_area_ratio
        land_fuel_area_ratio[shape] = land.fuel_area_ratio

        land_shape_index_biod[shape] = land.landscape_shape_index()
        land_shape_index_fuel[shape] = land.landscape_shape_index(var="fuel")

        land_diameter_biod[shape] = land.diameter()
        land_diameter_fuel[shape] = land.diameter(var='fuel')

        land_connectivity_correlation_biod[shape] = land.connectivity_correlation()
        land_connectivity_correlation_fuel[shape] = land.connectivity_correlation(var="fuel")

        land_IIC_nc_biod[shape] = land.IIC_with_nc()
        land_IIC_nc_fuel[shape] = land.IIC_with_nc(var="fuel")

        land_IIC_biod[shape] = land.IIC_without_nc()
        land_IIC_fuel[shape] = land.IIC_without_nc(var="fuel")

        land_CPL_nc_biod[shape] = land.CPL()
        land_CPL_nc_fuel[shape] = land.CPL(var='fuel')

        land_CPL_biod[shape] = land.CPL_without_nc()
        land_CPL_fuel[shape] = land.CPL_without_nc(var="fuel")

        components_number_biod[shape] = land.components_biod_tot[0]
        components_number_fuel[shape] = land.components_fuel_tot[0]

        components_area_biod_max[shape] = max(land.components_biod_tot[1])
        components_area_fuel_max[shape] = max(land.components_fuel_tot[1])

        components_shape_index_biod[shape] = land.components_shape_index()
        components_shape_index_fuel[shape] = land.components_shape_index("fuel")
        # print(values_unique.index(shape)/len(values_unique))

    checker["nodes_biod"] = [nodes_biod.get(n, n) for n in values]
    checker["nodes_fuel"] = [nodes_fuel.get(n, n) for n in values]
    checker["nodes_zero"] = [16 - checker.nodes_biod.loc[x] for x in range(len(checker))]
    checker["score_biod"] = [biod_score.get(n, n) for n in values]
    checker["score_fuel"] = [fuel_score.get(n, n) for n in values]

    checker["land_perimeter_biod"] = [land_perimeter_biod.get(n, n) for n in values]
    checker['land_perimeter_fuel'] = [land_perimeter_fuel.get(n, n) for n in values]
    checker["land_shape_index_biod"] = [land_shape_index_biod.get(n, n) for n in values]
    checker["land_shape_index_fuel"] = [land_shape_index_fuel.get(n, n) for n in values]

    checker["land_diameter_biod"] = [land_diameter_biod.get(n, n) for n in values]
    checker["land_diameter_fuel"] = [land_diameter_fuel.get(n, n) for n in values]

    checker["land_connectivity_correlation_biod"] = [land_connectivity_correlation_biod.get(n, n) for n in values]
    checker["land_connectivity_correlation_fuel"] = [land_connectivity_correlation_fuel.get(n, n) for n in values]

    checker["land_IIC_nc_biod"] = [land_IIC_nc_biod.get(n, n) for n in values]
    checker["land_IIC_nc_fuel"] = [land_IIC_nc_fuel.get(n, n) for n in values]

    checker["land_IIC_biod"] = [land_IIC_biod.get(n, n) for n in values]
    checker["land_IIC_fuel"] = [land_IIC_fuel.get(n, n) for n in values]

    checker["land_CPL_nc_biod"] = [land_CPL_nc_biod.get(n, n) for n in values]
    checker["land_CPL_nc_fuel"] = [land_CPL_nc_fuel.get(n, n) for n in values]

    checker["land_CPL_biod"] = [land_CPL_biod.get(n, n) for n in values]
    checker["land_CPL_fuel"] = [land_CPL_fuel.get(n, n) for n in values]

    checker["components_number_biod"] = [components_number_biod.get(n, n) for n in values]
    checker["components_number_fuel"] = [components_number_fuel.get(n, n) for n in values]

    checker["components_area_max_biod"] = [components_area_biod_max.get(n, n) for n in values]
    checker["components_area_max_fuel"] = [components_area_fuel_max.get(n, n) for n in values]

    checker["components_shape_index_biod"] = [components_shape_index_biod.get(n, n) for n in values]
    checker["components_shape_index_fuel"] = [components_shape_index_fuel.get(n, n) for n in values]

    convergence_time = []
    convergence_period = []
    indices = list(range(0, len(checker) + 1, 20))
    indices.remove(0)
    for x in indices:
        succession = checker.value.loc[x - 20:x - 1]
        result = functions_project.convergence(succession)
        convergence_time_s = result[0]
        convergence_period_s = result[1]
        convergence_time.extend([convergence_time_s] * 20)
        convergence_period.extend([convergence_period_s] * 20)
    del result

    checker['convergence_time'] = convergence_time
    checker['convergence_period'] = convergence_period
    checker['index_biod'] = [12] * len(checker)
    checker['initial_nb_1'] = [4] * len(checker)
    checker['initial_nb_2'] = [5] * len(checker)

    checker.to_csv('C:/Users/jean/Desktop/essai_stats.csv')
    print(time.time() - debut)

    times.append(time.time() - debut)
    timess+=1
