import numpy as np
import pandas as pd

from metpy.gridding.triangles import *
from metpy.gridding import points, triangles, polygons, gridding_functions, interpolation

from scipy.spatial import cKDTree, Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d, ConvexHull, qhull

import csv

def nn_point(xp, yp, grid_loc, tri, neighbors, triangle_info):
    r"""Generate a natural neighbor interpolation of the observations to the given point.

    This uses the Liang and Hale approach [Liang2010]_. The interpolation will fail if
    the grid point has no natural neighbors.

    Parameters
    ----------
    xp: (N, ) ndarray
        x-coordinates of observations
    yp: (N, ) ndarray
        y-coordinates of observations
    grid_loc: (float, float)
        Coordinates of the grid point at which to calculate the
        interpolation.
    tri: object
        Delaunay triangulation of the observations.
    neighbors: (N, ) ndarray
        Simplex codes of the grid point's natural neighbors. The codes
        will correspond to codes in the triangulation.
    triangle_info: dictionary
        Pre-calculated triangle attributes for quick look ups. Requires
        items 'cc' (circumcenters) and 'r' (radii) to be associated with
        each simplex code key from the delaunay triangulation.

    Returns
    -------
    value: float
       Interpolated value for the grid location

    """
    edges = triangles.find_local_boundary(tri, neighbors)
    # print 'edges:', edges
    edge_vertices = [segment[0] for segment in polygons.order_edges(edges)]
    num_vertices = len(edge_vertices)

    # edge_vertices are the indices of observation pts
    # print 'edge_vertices:', edge_vertices
    
    p1 = edge_vertices[0]
    p2 = edge_vertices[1]
#     print "p1:", p1
#     print "p2:", p2
    c1 = triangles.circumcenter(grid_loc, tri.points[p1], tri.points[p2])
    # print 'c1:', c1
    polygon = [c1]
#     print 'c1:',polygon
    area_list = []
    non_neighbor_weights = {}
    neighbor_weights = {}
    
    idx = range(len(xp))
    non_neighbors = [i for i in idx if i not in edge_vertices]
    for pt in non_neighbors:
        non_neighbor_weights['station_' + str(pt+1)] = 0

    total_area = 0.0
    for i in range(num_vertices):

        p3 = edge_vertices[(i + 2) % num_vertices]
        # print "p3:", p3
        try:
            # obtain index of neighbors
            pts_idx = edge_vertices.index(p3) - 1
            neighbor_idx = edge_vertices[pts_idx]
            
            c2 = triangles.circumcenter(grid_loc, tri.points[p3], tri.points[p2])
            polygon.append(c2)

            for check_tri in neighbors:
                if p2 in tri.simplices[check_tri]:
                    polygon.append(triangle_info[check_tri]['cc'])
            # print 'inbetween polygon:', polygon
            pts = [np.round(polygon[i],decimals=2) for i in ConvexHull(polygon).vertices]
            # print 'pts are:',pts

            cur_area = polygons.area(pts)
            neighbor_weights['station_' + str(neighbor_idx+1)] = cur_area
            
            # print 'cur_area___:', cur_area
            total_area += cur_area
            area_list.append(cur_area)

        except (ZeroDivisionError, qhull.QhullError) as e:
            message = ('Error during processing of a grid. '
                       'Interpolation will continue but be mindful '
                       'of errors in output. ') + str(e)

            log.warning(message)
            return np.nan

        polygon = [c2]
        # print 'polygon:', polygon
        p2 = p3
        # print 'new p2:', p2
    if total_area == sum(area_list):
        print 'True!'
        neighbor_weights.update((k, v/total_area) for k,v in neighbor_weights.items())
        total_weights = dict(neighbor_weights,**non_neighbor_weights)
        
        return total_weights
    else:
        print 'The area is not enclosed!'
        


#### RUN ####
# Input data from excel
cells    = pd.read_excel('cells_coordinates.xlsx', sheetname='lxly')
stations = pd.read_excel('wind_coordiantes.xlsx', sheetname='wind_coordiantes')

# Transform the dataframe into numpy arrays, only x, y coordiantes are needed
cells_coordinates    = cells[cells.columns[-2:]].get_values()
stations_coordinates = stations[stations.columns[-2:]].get_values()

xp  = stations_coordinates[:,0]
yp  = stations_coordinates[:,1]
pts = stations_coordinates

tri = Delaunay(pts)
# delaunay_plot_2d(tri)

    
sim_gridx = list(cells_coordinates[:, 0])
sim_gridy = list(cells_coordinates[:, 1])


members, tri_info = find_natural_neighbors(tri, list(zip(sim_gridx, sim_gridy)))

list_of_weights_dict = []

for i in range(len(sim_gridx)):
    edges = triangles.find_local_boundary(tri, members[i])
    if not edges:
    	i += 1
    	print i
        print edges
        weights = {}
        station_names = ['station_' + str(i+1) for i in range(len(pts))]
        for s_n in station_names:
        	weights[s_n] = 0 
        list_of_weights_dict.append(weights)
    else:
        weights = nn_point(xp, yp, [sim_gridx[i], sim_gridy[i]], tri, members[i], tri_info)
        list_of_weights_dict.append(weights)

# save list of dicts into csv
with open('nn_weights_15stations.csv', 'wb') as f: 
    station_names = ['station_' + str(i+1) for i in range(len(pts))]
    w = csv.DictWriter(f, list_of_weights_dict[0].keys())
    w.writeheader()
    w.writerows(list_of_weights_dict)