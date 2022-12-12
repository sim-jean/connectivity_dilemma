import ast
import os

import current_scripts.functions_project as functions_project
import current_scripts.outputs as results

exec(open('current_scripts/functions_project.py').read())

data = pd.read_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_2/successions_matched/land4_budget_2_9_5_cut_9_biodiv_index_12.csv").drop(columns=['Unnamed: 0'])

data               = data.iloc[:,0:20]
names              = list(range(20))
data.columns       = names

data["index"]      = list(range(len(data)))
data['index_biod'] = [12]*len(data)
#could be better to do in the end, with initial numbers as well.

checker            = data.melt(id_vars=["index","index_biod"])
checker["variable"] = checker['variable'].astype(int)
checker            = checker.sort_values(["index","variable"])
checker            = checker.reset_index(drop=True)

values             = list(checker.value)
values_unique      = list(set(values))

#region Setting up of data storage dictionnaries
nodes_biod         = {}
nodes_fuel         = {}

biod_score         = {}
fuel_score         = {}

land_perimeter_biod = {}
land_perimeter_fuel = {}

land_shape_index_biod = {}
land_shape_index_fuel = {}

land_diameter_biod = {}
land_diameter_fuel = {}

land_connectivity_correlation_biod = {}
land_connectivity_correlation_fuel = {}

land_IIC_nc_biod  = {}
land_IIC_nc_fuel  = {}

land_IIC_biod     = {}
land_IIC_fuel     = {}

land_CPL_nc_biod  = {}
land_CPL_nc_fuel  = {}

land_CPL_biod     = {}
land_CPL_fuel     = {}

components_number_biod = {}
components_number_fuel = {}

components_area_biod_max = {}
components_area_fuel_max = {}

components_shape_index_biod = {}
components_shape_index_fuel = {}
#endregion

#region Compute statistics
debut = time.time()
for shape in values_unique:
    land = Land(np.array(ast.literal_eval(shape)))

    nodes_biod[shape]    = land.node_biod
    nodes_fuel[shape]    = land.node_fuel

    biod_score[shape]    = land.biod
    fuel_score[shape]    = land.fuel

    land_perimeter_biod[shape] = land.perimeter()
    land_perimeter_fuel[shape] = land.perimeter(var="fuel")

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

    components_number_biod[shape] = land.components()
    components_number_fuel[shape] = land.components("fuel")

    components_area_biod_max[shape] = land.components_area_max()
    components_area_fuel_max[shape] = land.components_area_max("fuel")

    components_shape_index_biod[shape] = land.components_shape_index()
    components_shape_index_fuel[shape] = land.components_shape_index("fuel")
    #print(values_unique.index(shape)/len(values_unique))

checker["nodes_biod"] = [nodes_biod.get(n, n) for n in values]
checker["nodes_fuel"] = [nodes_fuel.get(n, n) for n in values]
checker["nodes_zero"] = [16 - checker.nodes_biod.loc[x] for x in len(checker)]
checker["score_biod"] = [biod_score.get(n, n) for n in values]
checker["score_fuel"] = [fuel_score.get(n, n) for n in values]


checker["land_perimeter_biod"]   = [land_perimeter_biod.get(n, n) for n in values]
checker['land_perimeter_fuel']   = [land_perimeter_fuel.get(n, n) for n in values]
checker["land_shape_index_biod"] = [land_shape_index_biod.get(n, n) for n in values]
checker["land_shape_index_fuel"] = [land_shape_index_fuel.get(n, n) for n in values]

checker["land_diameter_biod"]    = [land_diameter_biod.get(n, n) for n in values]
checker["land_diameter_fuel"]    = [land_diameter_fuel.get(n, n) for n in values]

checker["land_connectivity_correlation_biod"] = [land_connectivity_correlation_biod.get(n, n) for n in values]
checker["land_connectivity_correlation_fuel"] = [land_connectivity_correlation_fuel.get(n, n) for n in values]

checker["land_IIC_nc_biod"]       = [land_IIC_nc_biod.get(n, n) for n in values]
checker["land_IIC_nc_fuel"]       = [land_IIC_nc_fuel.get(n, n) for n in values]

checker["land_IIC_biod"]          = [land_IIC_biod.get(n, n) for n in values]
checker["land_IIC_fuel"]          = [land_IIC_fuel.get(n, n) for n in values]

checker["land_CPL_nc_biod"]       = [land_CPL_nc_biod.get(n, n) for n in values]
checker["land_CPL_nc_fuel"]       = [land_CPL_nc_fuel.get(n, n) for n in values]

checker["land_CPL_biod"]          = [land_CPL_biod.get(n, n) for n in values]
checker["land_CPL_fuel"]          = [land_CPL_fuel.get(n, n) for n in values]

checker["components_number_biod"] = [components_number_biod.get(n, n) for n in values]
checker["components_number_fuel"] = [components_number_fuel.get(n, n) for n in values]

checker["components_area_max_biod"] = [components_area_biod_max.get(n, n) for n in values]
checker["components_area_max_fuel"] = [components_area_fuel_max.get(n, n) for n in values]

checker["components_shape_index_biod"] = [components_shape_index_biod.get(n, n) for n in values]
checker["components_shape_index_fuel"] = [components_shape_index_fuel.get(n, n) for n in values]
print(time.time()-debut)

checker.to_csv('C:/Users/jean/Desktop/essai_stats.csv')
print(time.time()-debut)

