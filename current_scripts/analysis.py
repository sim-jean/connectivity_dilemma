import ast
import os

import numpy as np

import current_scripts.functions_project as functions_project
exec(open('current_scripts/functions_project.py').read())

import pandas as pd
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_2/successions_matched/land4_budget_2_10_0_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])

successions = data.iloc[550,0:15]

lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]

for x in lands:
    x.view()

new_dat = {}
new_dat['land']= [lands[x].shape for x in range(len(lands))]
new_dat['perim_biodiv']=[perimeter_biod(lands[x].mat) for x in range(len(lands))]
new_dat['biod_indic']= [lands[x].biod for x in range(len(lands))]
new_dat['fuel_indic']= [lands[x].fuel for x in range(len(lands))]
new_dat['perim_fire']=[perimeter_fire(lands[x].mat) for x in range(len(lands))]
new_dat['charact_path_biod']=[lands[x].characteristic_path_length_biod() for x in range(len(lands))]
new_dat['charact_path_fire']=[lands[x].characteristic_path_length_fire() for x in range(len(lands))]
new_dat['diameter_biod']=[lands[x].diameter_biod() for x in range(len(lands))]
new_dat['diameter_fire']=[lands[x].diameter_fire() for x in range(len(lands))]
new_dat['IIC_biod']=[lands[x].IIC_biod() for x in range(len(lands))]
new_dat['IIC_fuel']=[lands[x].IIC_fuel() for x in range(len(lands))]
data_eval = pd.DataFrame(new_dat)

# New metric that would be interesting is to verify whether an observed path of convergence is common. Like the checkered path, or the the parallel path
# First, evaluate when there is a convergence, or the length of a cycle.
# Find the first t after t_initial such that :
# nzero = budget?
#

data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_2/successions_matched/land4_budget_2_14_7_cut_0_biodiv_index_2.csv").drop(columns=['Unnamed: 0'])

successions = data.iloc[550,0:15]

lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]

for x in lands:
    x.view()

new_dat = {}
new_dat['land']= [lands[x].shape for x in range(len(lands))]
new_dat['perim_biodiv']=[perimeter_biod(lands[x].mat) for x in range(len(lands))]
new_dat['biod_indic']= [lands[x].biod for x in range(len(lands))]
new_dat['fuel_indic']= [lands[x].fuel for x in range(len(lands))]
new_dat['perim_fire']=[perimeter_fire(lands[x].mat) for x in range(len(lands))]
new_dat['charact_path_biod']=[lands[x].characteristic_path_length_biod() for x in range(len(lands))]
new_dat['charact_path_fire']=[lands[x].characteristic_path_length_fire() for x in range(len(lands))]
new_dat['diameter_biod']=[lands[x].diameter_biod() for x in range(len(lands))]
new_dat['diameter_fire']=[lands[x].diameter_fire() for x in range(len(lands))]
new_dat['IIC_biod']=[lands[x].IIC_biod() for x in range(len(lands))]
new_dat['IIC_fuel']=[lands[x].IIC_fuel() for x in range(len(lands))]
data_eval = pd.DataFrame(new_dat)

# What would be interesting is to compile those indicators and check how they respond to our lot of variables
# such as budget, or biodiversity value, and run analysis on the whole dataset. So the work is then to deal with those indicators and find the relevant ones.

#region check what's up with one
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/successions_matched/land4_budget_1_10_4_cut_2_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])

successions = data.iloc[550,0:15]

lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]

for x in lands:
    x.view()

new_dat = {}
new_dat['land']= [lands[x].shape for x in range(len(lands))]
new_dat['perim_biodiv']=[perimeter_biod(lands[x].mat) for x in range(len(lands))]
new_dat['biod_indic']= [lands[x].biod for x in range(len(lands))]
new_dat['fuel_indic']= [lands[x].fuel for x in range(len(lands))]
new_dat['perim_fire']=[perimeter_fire(lands[x].mat) for x in range(len(lands))]
new_dat['charact_path_biod']=[lands[x].characteristic_path_length_biod() for x in range(len(lands))]
new_dat['charact_path_fire']=[lands[x].characteristic_path_length_fire() for x in range(len(lands))]
new_dat['diameter_biod']=[lands[x].diameter_biod() for x in range(len(lands))]
new_dat['diameter_fire']=[lands[x].diameter_fire() for x in range(len(lands))]
new_dat['IIC_biod']=[lands[x].IIC_biod() for x in range(len(lands))]
new_dat['IIC_fuel']=[lands[x].IIC_fuel() for x in range(len(lands))]
data_eval = pd.DataFrame(new_dat)
#Need to store in the original nb2 nonz and zero.

data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/successions_matched/NEWland4_budget_1_10_4_cut_2_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])

successions = data.iloc[550,0:15]

lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]

for x in lands:
    x.view()


# Observation with 3
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_3/successions_matched/land4_budget_3_8_3_cut_0_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])

successions = data.iloc[550,0:15]

lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]

for x in lands:
    x.view()


data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/successions_matched/land4_budget_1_15_10_biodiv_index_2.csv").drop(columns=['Unnamed: 0'])
successions = data.iloc[0,0:15]
lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]

for x in lands:
    x.view()

#converges all right but does not match the biodiv index nor the initial conditions.

data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/land4_budget_1_2_2.csv")
#endregion
#region Where are we on the data?
lengths = []
for budget in [1,2,3,4]:
    files = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/successions_matched/"
    files = os.listdir(files)
    lengths.append(len(files)/(605*9))


#endregion


data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_2/successions_matched/land4_budget_2_9_5_cut_9_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])
successions = data.iloc[0,1:20]

lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]



horizon = list(range(len(successions)))
horizon.remove(0)
convergence_time = []
for t in horizon:
    if len(lands[t].equivalent().intersection(set([tuple(lands[t-1].shape)])))>0:
        convergence_time.append(t)
        break

for t in range(convergence_time[0]+1,horizon[-1],1):
    if tuple(lands[t].shape)==tuple(lands[convergence_time[0]].shape):
        cycle=t-convergence_time[0]
        break

#test

new_dat = {}
new_dat['land']= [lands[x].shape for x in range(len(lands))]
new_dat['perim_biodiv']=[perimeter_biod(lands[x].mat) for x in range(len(lands))]
new_dat['biod_indic']= [lands[x].biod for x in range(len(lands))]
new_dat['fuel_indic']= [lands[x].fuel for x in range(len(lands))]
new_dat['perim_fire']=[perimeter_fire(lands[x].mat) for x in range(len(lands))]
new_dat['charact_path_biod']=[lands[x].characteristic_path_length_biod() for x in range(len(lands))]
new_dat['charact_path_fire']=[lands[x].characteristic_path_length_fire() for x in range(len(lands))]
new_dat['diameter_biod']=[lands[x].diameter_biod() for x in range(len(lands))]
new_dat['diameter_fire']=[lands[x].diameter_fire() for x in range(len(lands))]
new_dat['IIC_biod']=[lands[x].IIC_biod() for x in range(len(lands))]
new_dat['IIC_fuel']=[lands[x].IIC_fuel() for x in range(len(lands))]
data_eval = pd.DataFrame(new_dat)

# Need to prepare pipeline for statistics implementation
# Do we need to keep the time series structure? Do we want to get at a panel with various levels of biodiv constraint?
# land - biodiv - t
# land - biodiv - t+1
# land - biodiv - t+2
# ...
# land - biodiv2 - var1 t  - var2 t
# land - biodiv2 - var1 t+1 - var2 t+1
debut = time.time()

data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_2/successions_matched/land4_budget_2_9_5_cut_9_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])

data                = data.iloc[:,0:20]
names               = list(range(20))
data.columns        = names

data["index"]       = list(range(len(data)))
data['index_biod']  = [12]*len(data)
#could be better to do in the end, with initial numbers as well.

checker             = data.melt(id_vars=["index","index_biod"])
checker["variable"] = checker['variable'].astype(int)
checker['land']     = checker.value.apply(lambda x: Land(shape=np.array(ast.literal_eval(x))))
checker             = checker.sort_values(["index","variable"])
checker             = checker.reset_index(drop=True)

checker['biodiv']   = checker.land.apply(lambda x : x.biod)
checker['fuel']     = checker.land.apply(lambda x : x.fuel)
checker['IIC_biod'] = checker.land.apply(lambda x: x.IIC_without_nc(var="biod"))
print(time.time()-debut)


a= checker.groupby("index").apply(lambda x: x.sort_values('variable'))
b= checker.sort_values("index").groupby("variable").head(10)

debut = time.time()
essai = list(checker.value)
essai_unique = list(set(essai))
#initiate variable
IIC_biod = []
IIC_fuel = []

diam_biod = []
diam_fuel = []
for x in essai_unique:
    a = Land(shape=np.array(ast.literal_eval(x)))
    IIC_biod.append(a.IIC_without_nc)
    IIC_fuel.append(a.IIC_without_nc(var="fuel"))
    diam_biod.append(a.diameter)
    diam_fuel.append(a.diameter(var="fuel"))
    print(str(essai_unique.index(x)/len(essai_unique)))
print(time.time()-debut)

# 1 min 30 au moins pour chaque dataset. Disons 2 avec les autres métriques.
total_per_budg = 605*9
total = total_per_budg*4
time_per_run = 5
total_length = total*time_per_run
total_duration_mp = total_length/16
total_duration_mp/(60*24)

# Variables
IIC_biod = {}
IIC_fuel = {}

diam_biod = {}
diam_fuel = {}

CPL_biod = {}
CPL_fuel = {}

connectivity_correlation_biod = {}
connectivity_correlation_fuel = {}

perimeter_biod = {}
perimeter_fuel = {}

perimeter_to_area_biod = {}
perimeter_to_area_fuel = {}

nb_component_biod = {}
nb_component_fuel = {}

mc_perimeter_biod = {}
mc_perimeter_fuel = {}

mc_perimeter_to_area_biod = {}
mc_perimeter_to_area_fuel = {}

mc_diam_biod = {}
mc_diam_fuel = {}

habitat ={}
fuel = {}

biod_index = {}
fuel_index = {}

components_area_biod = {}
components_area_fuel = {}
debut = time.time()
essai = list(checker.value)
essai_unique = list(set(essai))
#initiate variable
for x in essai_unique:
    a = Land(shape=np.array(ast.literal_eval(x)))
    diam_fuel[x]=a.diameter(var="fuel")
    diam_biod[x]=a.diameter(var="biod")

    IIC_biod[x] = a.IIC_without_nc(var="biod")
    IIC_fuel[x] = a.IIC_without_nc(var="fuel")

    CPL_biod[x] = a.characteristic_path_length(var='biod')
    CPL_fuel[x] = a.characteristic_path_length(var='fuel')

    connectivity_correlation_biod[x] = a.connectivity_correlation(var='biod')
    connectivity_correlation_fuel[x] = a.connectivity_correlation(var='fuel')

    nb_component_biod[x] = a.components(var="biod")[0]
    nb_component_fuel[x] = a.components(var="fuel")[0]

    perimeter_biod[x] = a.perimeter(var="biod")
    perimeter_fuel[x] = a.perimeter(var="fuel")

    habitat[x] = sum(a.shape>=1)
    fuel[x] = sum(a.shape==2)

    biod_index[x]=a.biod
    fuel_index[x]=a.fuel

    components_area_biod[x]=a.components_area(var="biod")
    components_area_fuel[x]=a.components_area(var="fuel")


    #print(str(essai_unique.index(x)/len(essai_unique)))
print(time.time()-debut)



# verify if data is coherent, and maybe have to use hyperindexing
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_2/successions_matched/land4_budget_2_9_5_cut_3_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])
budg2 = data.iloc[:,0]
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/successions_matched/land4_budget_1_9_5_cut_3_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])
budg1 = data.iloc[:,0]
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_3/successions_matched/land4_budget_3_9_5_cut_3_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])
budg3 = data.iloc[:,0]

##"
budg1

succession_dataframe(data_path,12)
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_1/successions_matched/land4_budget_1_9_5_cut_2_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])
budg1 = data.iloc[:,0]
successions = data.iloc[0,0:20]
successions = [Land(np.array(ast.literal_eval(x))) for x in successions]

successions[0].view()
successions[1].view()
successions[2].view()
successions[3].view()
successions[4].view()
successions[5].view()
successions[6].view()
successions[7].view()
successions[8].view()
successions[9].view()
successions[10].view()

data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_3/successions_matched/land4_budget_3_10_4_cut_4_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])
budg1 = data.iloc[:,0]
successions = data.iloc[0,0:20]
successions = [Land(np.array(ast.literal_eval(x))) for x in successions]

successions[0].view()
successions[1].view()
successions[2].view()
successions[3].view()
successions[4].view()
successions[5].view()
successions[6].view()
successions[7].view()
successions[8].view()
successions[9].view()
successions[10].view()


land
land2 = to_matrix4(land)
land2 = [land2[3],land2[2],land2[1],land2[0]]
land2 = tuple(list(itertools.chain(*land2)))
land2 = np.rot90(np.rot90(np.rot90(to_matrix4(land2))))
land2 = tuple(list(itertools.chain(*land2)))

land3 = np.rot90(np.rot90(np.rot90(land2,axes=(1,0)),axes=(1,0)),axes=(1,0))
land3 = [land3[3],land3[2],land3[1],land3[0]]
land3 = tuple(list(itertools.chain(*land3)))
land3 = to_matrix4(land3)
land3==to_matrix4(land)


#shapes
shape4 = np.array((1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1))
shape5 = np.array((0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0))

land4 = Land(shape4)
land4.view()
land4.biod
land5 = Land(shape5)
land5.view()
land5.biod
land6 = Land(shape= np.array((0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0)))
land6.view()
land6.biod
land7 = Land(shape = np.array(ast.literal_eval(data.iloc[458,0])))
land7.view()
land7.fuel
land7.biod

land4b = Land(np.array((1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1)))
land4b.view()
land4b.biod
land4.biod
land4c =  Land(np.array((0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1)))
land4c.view()


shape = np.repeat(2,16)
land = Land(shape)
pot_fire_budget = [x for x in pot_fire if sum(sum(to_matrix(x,size)*costs))<=5]
land.fuel_dynamics(pot_fire_budget)
shape2 = land.succession([2])[0][2]
land2 = Land(shape2)
land2.view()

land2.fuel_dynamics(pot_fire_budget)
land2.succession([2])
land2.view_succession([2])
shape3 = land2.succession([2])[0][2]
land3 = Land(shape3)
land3.view()
land3.fuel_dynamics(pot_fire_budget)
land3.succession([2])
land3.view_succession([2])
shape4 = land3.succession([2])[0][2]
land4 = Land(shape4)
land4.fuel_dynamics(pot_fire_budget)
land4.succession([2])
land4.view_succession([2])
shape5 = land4.succession([2])[0][2]
land5 = Land(shape5)
land5.fuel_dynamics(pot_fire_budget)
land5.succession([2])
land5.view_succession([2])

# try for a specific cut
sample = "land4_budget_1_10_5_cut_15"
data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_1/successions_matched/"+sample+"_biodiv_index_82.csv")
shape = np.array(ast.literal_eval(data.iloc[0,19]))
land = Land(shape)
land.fuel
tryer={"value":[],"biod":[],"ones":[],"twos":[]}

files = os.listdir("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_3/successions_matched/")
for file in files:
    data = pd.read_csv(
    "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_3/successions_matched/" +file)
    z = re.findall(r'\d+', file)
    data_append = list(set(data.iloc[:,41]))
    tryer["value"].extend(data_append)
    tryer["biod"].extend([z[-1]]*len(data_append))
    tryer["ones"].extend([int(z[2])-int(z[3])]*len(data_append))
    tryer["twos"].extend([z[3]]*len(data_append))
    print(files.index(file)/len(files))

tester = pd.DataFrame(tryer)
tester["og_data"]=tester[["ones","twos"]].apply(tuple,axis=1)

def dfScatter(df, xcol='biod', ycol='value', catcol='og_data'):
    fig, ax = plt.subplots()
    categories = np.unique(df[catcol])
    colors = np.linspace(0, 1, len(categories))
    colordict = dict(zip(categories, colors))

    df["Color"] = df[catcol].apply(lambda x: colordict[x])
    ax.scatter(df[xcol], df[ycol], c=df.Color)
    return fig

fig = dfScatter(tester)
fig.savefig("fig1.png")