#################################
# Check if successions are well paired
import ast
import os
import random

import pandas as pd

budget=1
path = 'C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_'+str(budget)+'/successions_matched/'
files = os.listdir(path)

# budget = 1
data_check = pd.read_csv(path+files[2])

successions = data_check.iloc[1,1:11]
lands = []
for x in successions:
    lands.append(Land(shape=np.array(ast.literal_eval(x))))

for x in lands:
    x.view()

# Does not work with budget 1....
index = random.randint(0,len(files))
data_check = pd.read_csv(path+files[index])
data_check2 = pd.read_csv(path+'NEWland4_budget_1_10_5_cut_23_biodiv_index_2.csv')
data_check = pd.read_csv(path+"land4_budget_1_10_5_cut_23_biodiv_index_72.csv")

successions = data_check.iloc[0,1:11]
successions2 = data_check2.iloc[0,1:11]
lands = []
for x in successions:
    lands.append(Land(shape=np.array(ast.literal_eval(x))))
lands[0].view()
lands[1].view()
lands[2].view()
lands[3].view()
lands[4].view()
lands[5].view()

lands2 = []
successions3 = [ast.literal_eval(x) for x in successions2]
for x in successions3:
    lands2.append(Land(shape=np.array(x)))
lands2[0].view()
lands2[1].view()
lands2[2].view()
lands2[3].view()
lands2[4].view()
lands2[5].view()
lands2[6].view()
index_biod=2

def succession_dataframe(data_path, biodiv_index,time_param='Yes',budget=3,follow_up="No"):

    debut = time.time()
    #region Load the current data

    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(16 - budget, 17))]
    datas = {}
    datas[data_path]=rows_data
    for x in files:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/" +x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/"
    files = os.listdir(path)
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[0]) in list(range(16 - budget, 17))]
    datas_equiv = {}
    for x in files:
        # region Load data in list format
        file = open(path + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[x] = rows_dict
        #endregion
    # endregion
    #region Add current data to dict repositories

    match_equiv = re.findall(r'\d+', str(data_path))
    if len(match_equiv)>4:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+'.csv')]=rows_dict
    else:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+'.csv')]=rows_dict

    #endregion

    rows_out = []

    for p in list(range(len(rows_data))):
    # input
        land = str(rows_data[p][0])
    # output
        succession_list = [land]
        transformation = ['None']
        equiv_key = [land]
        value = [np.nan]
    #Meta parameters
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
        step = 0
        while step <= 20:
            check = nonzero(np.array(ast.literal_eval(equiv_key[-1])))
            check2 = nb2(np.array(ast.literal_eval(equiv_key[-1])))
            #print(check, check2)
    # use relevant result data from loaded data :
            namer = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
            namer2 = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + "_cut_0.csv")

            if os.path.exists(path + namer):
                if namer in list(datas.keys()):
                    d2 = datas[namer]
                else:
                    file = open(path + "land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
                    csvreader = csv.reader(file)
                    header = next(csvreader)

                    rows_dict = []
                    for row in csvreader:
                        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
                        rows_dict.append(row_add)
                    datas["land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv"] = rows_dict
                    d2 = rows_dict

                for g in list(range(len(d2))):
                    if d2[g][0]==equiv_key[-1]:
                        succession_index=g

                name2 = "equiv_results_" + str(check) + "_" + str(check2) + ".csv"

            elif os.path.exists(path + namer2):
                if namer2 in list(datas.keys()):
                    to_load = [x for x in list(datas.keys()) if x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                else :
                    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(
                        budget) + "/"
                    files = os.listdir(path)
                    files.remove('keys_succession')
                    files.remove('successions_matched')
                    to_load = [x for x in files if
                               x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                    for y in to_load:
                        file = open(path + y)
                        csvreader = csv.reader(file)
                        header = next(csvreader)

                        rows_dict = []
                        for row in csvreader:
                            row_add = [row[0], row[int(biodiv_index / 2)], row[int(82 + biodiv_index / 2)]]
                            rows_dict.append(row_add)
                        datas[y] = rows_dict

                for x in to_load:
                    d2 = datas[x]
                    for g in list(range(len(d2))):
                        if d2[g][0]==equiv_key[-1]:
                            succession_index = g
                            indexes = re.findall(r'\d+', str(x))
                            keep = x
                            name2 = "equiv_results_" + str(check) + "_" + str(check2) + "_cut_" + str(
                                indexes[4]) + ".csv"

                d2=datas[keep]
            succession_key = d2[succession_index][1]
            value_column = int(82 + biodiv_index / 2)

            if name2 in list(datas_equiv.keys()):
                d_equiv=datas_equiv[name2]
            else:
                file = open(path + "keys_succession/"+name2)
                csvreader = csv.reader(file)
                header = next(csvreader)

                rows_dict = []
                for row in csvreader:
                    row_add = row[int(biodiv_index/2-1)]
                    rows_dict.append(row_add)
                d_equiv = rows_dict
                datas_equiv[name2] = d_equiv

            #region Computations
            candidate_succ = d_equiv[succession_index]
            mat_succession = to_matrix4(ast.literal_eval(succession_key))
            transfo = equi_landscape4_recover(mat_succession).index(ast.literal_eval(candidate_succ))
            #endregion
            #region  list updates
            value.append(d2[succession_index][2])
            succession_list.append(succession_key)
            equiv_key.append(candidate_succ)
            transformation.append(transfo)
            #endregion

            step += 1
        #region Final processing : recovery of optimal succession with chain transformation
        c = list(map(ast.literal_eval, succession_list))
        final_succession = []
        #for i in list(range(20)):
        for i in list(range(5)):
            e = transformation[0: i ]
            # deleted +1
            final_succession.append(recover_full2(c[i], transfo=e, index=i))

        final_succession = final_succession + value
        rows_out.append(final_succession)
        if follow_up == 'Yes':
            print(str(100*len(rows_out) / len(rows_data))+" %")
        #endregion

    data_path_saved = data_path.replace('.csv', '_biodiv_index_' + str(biodiv_index) + '.csv')
    data_return = pd.DataFrame(rows_out)
    data_return.to_csv(path+'successions_matched/'+data_path_saved)


    if time_param == 'Yes':
        return(data_path + " is done and took "+str(time.time() - debut)+' seconds' )
    else :
        return(data_path + " is done")




### Budget 2
budget=2

path = 'C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_'+str(budget)+'/successions_matched/'
files = os.listdir(path)

data_check = pd.read_csv(path+files[2])

successions = data_check.iloc[0,1:11]
lands = []
for x in successions:
    lands.append(Land(shape=np.array(ast.literal_eval(x))))

for x in lands:
    x.view()

lands[0].view()
lands[1].view()
lands[2].view()
lands[3].view()

index = random.randint(0,len(files))
data_check = pd.read_csv(path+files[index])

successions = data_check.iloc[0,1:11]
lands = []
for x in successions:
    lands.append(Land(shape=np.array(ast.literal_eval(x))))
lands[0].view()
lands[1].view()
lands[2].view()
lands[3].view()
lands[4].view()
lands[5].view()

#
files ="C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_4/"
files = os.listdir(files)
files2 = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_4/"
files2 = os.listdir(files2)
files2.remove("successions_matched")
files2.remove("keys_succession")
#files2.remove("updated_list.csv")

not_already_done=[]
for x in files2:
    for biod in list(range(2,100,step)):
        y=x.replace(".csv",("_biodiv_index_"+str(biod)))
        y+=".csv"
        not_already_done.append(y)

not_already_done_set = set(not_already_done)
files_set = set(files)

c = not_already_done_set-files_set

d= list(c)
len(d)
tryer = [x.replace("_biodiv_index_"," ") for x in d]
tryer = [x.replace(" ", "_csv ") for x in tryer]
tryer = [x.replace(".csv","") for x in tryer]
tryer = [x.replace("_csv",".csv") for x in tryer]
tryer = [x.split(" ") for x in tryer]

candidates = []
for x in tryer:
    bc = (x[0],int(x[1]))
    candidates.append(bc)

with open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/updated_list.csv", 'w', newline='\n') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(candidates)

file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/updated_list.csv", "r")
data = list(csv.reader(file, delimiter=","))
file.close()
files = [ast.literal_eval(x) for x in data[0]]



candidates_flat=[ast.literal_eval(x) for x in data[0]]
step = 10
print(budget)
__name__='__main__'
__spec__=None
#candidates = []
#for j in files:
#    candidates.append(list(zip([j] * len(list(range(2, 84, step))), list(range(2, 84, step)))))


#candidates_flat = [item for sublist in candidates for item in sublist]
if __name__=="__main__":
    start = time.time()
    with mp.Pool(processes=16) as pool:
        results = pool.starmap(functions_project.succession_dataframe,candidates_flat)
    pool.terminate()
    print("The process took " + str((time.time()-start)//60) + ' minutes and ' + str((time.time()-start)%60) + ' seconds')





#
source = (x,y)
candidates = np.array(([False, False, False, False,],
                      [False, False, False, False,],
                      [False, False, False, False,],
                      [False, False, False, False,]))


# Program to count islands in boolean 2D matrix
class Graph:

    def __init__(self, row, col, g):
        self.ROW = row
        self.COL = col
        self.graph = g

    # A function to check if a given cell
    # (row, col) can be included in DFS
    def isSafe(self, i, j, visited):
        # row number is in range, column number
        # is in range and value is 1
        # and not yet visited
        return (i >= 0 and i < self.ROW and
                j >= 0 and j < self.COL and
                not visited[i][j] and self.graph[i][j])

    # A utility function to do DFS for a 2D
    # boolean matrix. It only considers
    # the 8 neighbours as adjacent vertices

    def DFS(self, i, j, visited):

        # These arrays are used to get row and
        # column numbers of 8 neighbours
        # of a given cell
        rowNbr = [-1, -1, -1, 0, 0, 1, 1, 1]
        colNbr = [-1, 0, 1, -1, 1, -1, 0, 1]

        # Mark this cell as visited
        visited[i][j] = True

        # Recur for all connected neighbours
        for k in range(8):
            if self.isSafe(i + rowNbr[k], j + colNbr[k], visited):
                self.DFS(i + rowNbr[k], j + colNbr[k], visited)

    # The main function that returns
    # count of islands in a given boolean
    # 2D matrix

    def countIslands(self):
        # Make a bool array to mark visited cells.
        # Initially all cells are unvisited
        visited = [[False for j in range(self.COL)]for i in range(self.ROW)]

        # Initialize count as 0 and traverse
        # through the all cells of
        # given matrix
        count = 0
        area = [0]
        for i in range(self.ROW):
            for j in range(self.COL):
                # If a cell with value 1 is not visited yet,
                # then new island found
                if visited[i][j] == False and self.graph[i][j] == 1:
                    # Visit all cells in this island
                    # and increment island count
                    self.DFS(i, j, visited)
                    count += 1
                    area.append(sum([sum(x) for x in visited])-area[-1])
        area.pop(0)

        return count, area


graph = [[1, 1, 0, 0, 0],
        [0, 1, 0, 0, 1],
        [1, 0, 0, 1, 1],
        [0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1]]


row = len(graph)
col = len(graph[0])

g = Graph(row, col, graph)

print("Number of islands is:")
print(g.countIslands())

# This code is contributed by Neelam Yadav
def isValid(x, y, grid, visited):
    if ((x >= 0 and y >= 0) and
            (x < len(grid) and y < len(grid[0])) and
            (grid[x][y] != 0) and (visited[x][y] == False)):
        return True
first_nodes = [(0,0),(0,3)]
if len(first_nodes)==1:
    land.perimeter()
else:
    components = []
    perimeter_components = []
    for nodes in first_nodes:
        queue = []
        queue.append(nodes)
        visited = np.ones((4,4),dtype=bool)==False
        patches = []
        cols = [0,0,-1,1,-1,-1,1,1]
        rows = [-1,1,0,0,1,-1,1,-1]
        while len(queue)>0:
            node = queue.pop(0)

            for i,j in list(zip(rows,cols)):
                if isValid(node[0]+i,node[1]+j,land.mat,visited):
                    queue.append((node[0]+i,node[1]+j))
                    visited[node[0]+i,node[1]+j]=True
                    patches.append((node[0]+i,node[1]+j))
        components.append(patches)
        perimeter=0
        for patch in patches:
            perimeter+=(land.mat[patch[0],patch[1]]>=1)*land.number_empty_neighbors(patch[0],patch[1])

        perimeter_components.append(perimeter)

    perimeter_components
    components


### New class Land_light()
# Same inits

class Land:
    def __init__(self,shape):
        self.shape= shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = to_matrix(shape,size)
        self.biod = (shape>=1).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=1)))-sum(shape>=1)
        self.biod_relative = self.biod/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape))))-sum(np.repeat(1,len(self.shape))))
        self.fuel = (shape>=2).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=2)))-sum(shape>=2)
        self.fuel_relative = self.fuel/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape))))-sum(np.repeat(1,len(self.shape))))
        self.node_biod = sum(shape>0)
        self.node_fuel = sum(shape==2)
        self.zero = sum(shape==0)
        self.habitat_area_ratio = sum(self.shape>=1)/(self.size**2)
        self.fuel_area_ratio = sum(self.shape==2)/(self.size**2)
        self.test_algo_0 = self.test_algo(0)
        self.test_algo_1 = self.test_algo(1)
        self.test_algo_2 = self.test_algo(2)
    def view(self):
        ax = plt.axes()
        ax.set_title('Initial biodiv = '+str(self.biod))
        sns.heatmap(self.mat,ax = ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
        plt.show()
    def view_succession(self,index):
        if len(index)==1:
            ax = plt.axes()
            ax.set_title('Biodiv constraint = ' + str(index[0])+" with budget = "+str(budget))
            y = to_matrix(self.succession_land[index[0]], self.size)
            sns.heatmap(y, ax=ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
            plt.show()
        else:
            number_x = 0
            rows = len(index)//2
            for i in range(rows):
                for j in range(2):
                    ax = plt.subplot2grid((rows, 2), (i, j))
                    ax.title.set_text("Biodiv constraint =" + str(index[number_x]) + ' with budget ='+str(budget))
                    y = to_matrix(self.succession_land[index[number_x]], self.size)
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
            biod =np.mat(np.diagonal((self.fuel_dynamics >= 1).dot(connectivity_mat(self.size)).dot(np.transpose((self.fuel_dynamics >= 1))) - (self.fuel_dynamics >= 1).sum(
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
            print("Problem")
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
    def diagonal(self,direction,var="biod"):
        if var=="biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)

        if direction == 'NO-SE':
            return(paths[((0,0),(self.size-1,self.size-1))])
        elif direction=='SO-NE':
            return(paths[((size-1,0),(0,size-1))])
    # Description
    def number_empty_neighbors(self,i,j,var="biod"):
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

        return g.countIslands()[0]
    def components_area(self,var='biod'):
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)
        result = g.countIslands()[1]
        if len(result)==1:
            result = int(result[0])
        elif len(result)==0:
            result = float("nan")
        return result
    def components_area_max(self,var='biod'):
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        result = g.countIslands()[1]
        if len(result) == 1:
            result = int(result[0])
        elif len(result) == 0:
            result = float("nan")
        else:
            result = max(result)
        return result
    def components_perimeter(self, var= "biod",components_var=False):
        vars = var
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)

        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row,col,graph)
        first_nodes = g.countIslands()[2]

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
            if self.components()>1:
                area = max(self.components_area(var="biod"))
                candidate = self.components_area(var="biod").index(area)
                perimeter = self.components_perimeter(var="biod")[candidate]
            else:
                area = self.components_area(var="biod")
                perimeter = self.components_perimeter(var="biod")

        elif var=="fuel":
            if self.components("fuel")>1:
                area = max(self.components_area(var="fuel"))
                candidate = self.components_area(var="fuel").index(area)
                perimeter = self.components_perimeter(var="fuel")[candidate]
            else:
                area = self.components_area(var="fuel")
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
    def coord_zero(self):
        coord_zero=[]
        for i in range(self.size):
            for j in range(self.size):
                if self.mat[i,j]==0:
                    coord_zero.append((i,j))
        return coord_zero

def succession_dataframe(data_path, biodiv_index,time_param='Yes',budget=1,follow_up="No"):

    debut = time.time()
    #region Load the current data

    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    try:
        files.remove("updated_list.csv")
    except ValueError:
        pass

    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(16 - budget, 17))]
    datas = {}
    datas[data_path]=rows_data
    for x in files:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/" +x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/"
    files = os.listdir(path)
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[0]) in list(range(16 - budget, 17))]
    datas_equiv = {}
    for x in files:
        # region Load data in list format
        file = open(path + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[x] = rows_dict
        #endregion
    # endregion
    #region Add current data to dict repositories

    match_equiv = re.findall(r'\d+', str(data_path))
    if len(match_equiv)>4:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+'.csv')]=rows_dict
    else:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+'.csv')]=rows_dict

    #endregion

    rows_out = []

    for p in list(range(len(rows_data))):
    # input
        land = str(rows_data[p][0])
    # output
        succession_list = [land]
        transformation = ['None']
        equiv_key = [land]
        value = [np.nan]
    #Meta parameters
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
        step = 0
        while step <= 20:
            check = nonzero(np.array(ast.literal_eval(equiv_key[-1])))
            check2 = nb2(np.array(ast.literal_eval(equiv_key[-1])))
            #print(check, check2)
    # use relevant result data from loaded data :
            namer = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
            namer2 = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + "_cut_0.csv")

            if os.path.exists(path + namer):
                if namer in list(datas.keys()):
                    d2 = datas[namer]
                else:
                    file = open(path + "land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
                    csvreader = csv.reader(file)
                    header = next(csvreader)

                    rows_dict = []
                    for row in csvreader:
                        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
                        rows_dict.append(row_add)
                    datas["land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv"] = rows_dict
                    d2 = rows_dict
                #Find the right row that matches the equivalent landscape to original landscape (unicity)
                for g in list(range(len(d2))):
                    if d2[g][0]==equiv_key[-1]:
                        succession_index=g

                name2 = "equiv_results_" + str(check) + "_" + str(check2) + ".csv"

            elif os.path.exists(path + namer2):
                if namer2 in list(datas.keys()):
                    to_load = [x for x in list(datas.keys()) if x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                else :
                    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(
                        budget) + "/"
                    files = os.listdir(path)
                    files.remove('keys_succession')
                    files.remove('successions_matched')
                    to_load = [x for x in files if
                               x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                    for y in to_load:
                        file = open(path + y)
                        csvreader = csv.reader(file)
                        header = next(csvreader)

                        rows_dict = []
                        for row in csvreader:
                            row_add = [row[0], row[int(biodiv_index / 2)], row[int(82 + biodiv_index / 2)]]
                            rows_dict.append(row_add)
                        datas[y] = rows_dict

                for x in to_load:
                    d2 = datas[x]
                    for g in list(range(len(d2))):
                        if d2[g][0]==equiv_key[-1]:
                            succession_index = g
                            indexes = re.findall(r'\d+', str(x))
                            keep = x
                            name2 = "equiv_results_" + str(check) + "_" + str(check2) + "_cut_" + str(
                                indexes[4]) + ".csv"

                d2=datas[keep]

            # Get the succession to the equivalent/unique landscape.
            succession_key = d2[succession_index][1]
            value_column = int(82 + biodiv_index / 2)

            if name2 in list(datas_equiv.keys()):
                d_equiv=datas_equiv[name2]
            else:
                file = open(path + "keys_succession/"+name2)
                csvreader = csv.reader(file)
                header = next(csvreader)

                rows_dict = []
                for row in csvreader:
                    row_add = row[int(biodiv_index/2-1)]
                    rows_dict.append(row_add)
                d_equiv = rows_dict
                datas_equiv[name2] = d_equiv

            #region Computations
            #store the unique equivalent to the succession
            candidate_succ = d_equiv[succession_index]

            mat_succession = to_matrix4(ast.literal_eval(succession_key))
            # Based on the value of the succession of the unique landscape, find how this maps to the unique equivalent
            # to succession
            transfo = equi_landscape4_recover(mat_succession).index(ast.literal_eval(candidate_succ))
            #endregion
            #region  list updates
            value.append(d2[succession_index][2])
            succession_list.append(succession_key)
            equiv_key.append(candidate_succ)
            transformation.append(transfo)
            #endregion

            step += 1
        #region Final processing : recovery of optimal succession with chain transformation
        c = list(map(ast.literal_eval, succession_list))
        final_succession = []
        for i in list(range(20)):
            e = transformation[0: i+1 ]
            # deleted +1
            final_succession.append(recover_full2(c[i], transfo=e, index=i-1))
        final_succession[0]=ast.literal_eval(land)


        final_succession = final_succession + value
        rows_out.append(final_succession)
        if follow_up == 'Yes':
            print(str(100*len(rows_out) / len(rows_data))+" %")
        #endregion

    data_path_saved = data_path.replace('.csv', '_biodiv_index_' + str(biodiv_index) + '.csv')
    data_return = pd.DataFrame(rows_out)
    data_return.to_csv(path+'successions_matched/'+data_path_saved)


    if time_param == 'Yes':
        return(data_path + " is done and took "+str(time.time() - debut)+' seconds' )
        return(data_return)
    else :
        return(data_path + " is done")

data = pd.read_csv('C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_1/successions_matched/land3_budget_1_5_biodiv_index_2.csv')

successions = list(data.iloc[0,1:19])
successions2
def convergence(successions):
    """
    Returns the time of convergence and length of cycle if applicable
    :param successions: data series
    :return: a list of tuples (t_convergence, cycle_period)
    """
    lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]
    horizon = list(range(len(successions)))
    horizon.remove(0)
    convergence_time = -1
    for t in horizon:
        for i in range(1,t):
            if len(lands[t].equivalent().intersection(lands[t - i].equivalent())) > 0:
                break
        else:
            continue
        break
    convergence_time = t
    initial = i
    cycle = t-i

    for period in range(initial,horizon[-1]-cycle,cycle):
        if len(lands[period].equivalent().intersection(lands[period+cycle].equivalent()))==0:
            #print('Not a convergence')
            return float('nan'), float('nan')
            break

    else:
        #print('Convergence in ' + str(cycle))
        return convergence_time, cycle






lands[0].view()
lands[1].view()
lands[2].view()
lands[3].view()
lands[4].view()
lands[5].view()
lands[6].view()