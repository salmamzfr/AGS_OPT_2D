'''
Created on 13.04.2018, updated on 21.01.2020

@author: Salma Mozaffari, ETH Zurich

Comment: Helper Functions used in AGS_OPT_2D.py file
'''

"################################################ IMPORTS ###################################################"
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pickle
import numpy as np
from scipy.spatial import ConvexHull
from compas.datastructures import Mesh 
from compas.datastructures.network import Network
from compas.datastructures import network_find_crossings
from compas.datastructures import network_find_faces
from compas.geometry import Polygon
from compas.geometry import Polyline
from compas.geometry import is_point_in_polygon_xy
from compas.geometry import is_point_on_polyline
from compas.geometry import rotate_points
from compas.numerical.matrices import connectivity_matrix
from compas.utilities import geometric_key
from compas_plotters import MeshPlotter
from compas_plotters import NetworkPlotter

"##################################################################################################################"
"################################################ HELPER FUNCTIONS ################################################"
"##################################################################################################################"

"########################################## GENERAL FUNCTIONS #############################################"

def vector_create(pt_1_coor, pt_2_coor):
    """ 
    Returns numpy vector between 
    vertex 1 and vertex 2 coordinates
    """
    arr_1=np.array(pt_1_coor)
    arr_2=np.array(pt_2_coor)
    vec=arr_2-arr_1
    
    return vec

def unit_vector(vec):
    """
    returns the unit vector of vec array.
    """
    vec_norm=np.linalg.norm(vec)
    unit_vec=vec/vec_norm
    return unit_vec

def add_vectors(vec_1, vec_2):
    """
    adds two vectors: non-array results
    vec_1 & vec_2: XYZ components of the vectors.
    returns: list of resulting  vector
    """
    return [a+b for (a, b) in zip(vec_1, vec_2)]

def subtract_vectors(vec_1, vec_2):
    """
    subtract two vectors: non-array results
    vec_1 & vec_2: XYZ components of the vectors.
    returns: list of resulting  vector
    """
    return [a-b for (a, b) in zip(vec_1, vec_2)]

def middle_point(pt_1, pt_2):
    """
    finds the middle point of a line
    """
    sub=add_vectors(pt_1, pt_2)
    
    return (round(sub[0]/2.0, 4), round(sub[1]/2.0, 4), sub[2]/2.0)

def angle_between(vec_1, vec_2):
    """ 
    Returns the angle in radians between arrays vec_1 and vec_2
    """
    v1_u=unit_vector(vec_1)
    v2_u=unit_vector(vec_2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u),-1.0, 1.0))
    
def merge_two_dicts(x, y):
    """
    from Stackoverflow
    """
    z=x.copy()  # start with x's keys and values
    z.update(y)  # modifies z with y's keys and values & returns None
    return z

def distance_between_two_points(pt_1_coor, pt_2_coor):
    """
    returns the distance between two points 
    """
    vec=np.array(pt_1_coor)-np.array(pt_2_coor)
    vec_norm=np.linalg.norm(vec)
    return round(vec_norm, 5)

def edge_length(network, edg):
    """
    returns an edge length of a network
    """
    coor_1=network.vertex_coordinates(edg[0])
    coor_2=network.vertex_coordinates(edg[1])
    edg_len=distance_between_two_points(coor_1, coor_2)

    return edg_len
    
def line_parameters_xy(pt_1, pt_2):
    """
    Used in lines_intersection
    from:
    https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python
    """
    a=(pt_1[1]-pt_2[1])
    b=(pt_2[0]-pt_1[0])
    c=(pt_1[0]*pt_2[1]-pt_2[0]*pt_1[1])
    return a, b,-c

def lines_intersection_xy(line_1, line_2):
    """
    Cramer's Rule
    from:
    https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python
    """
    par_1=line_parameters_xy(line_1[0], line_1[1])
    par_2=line_parameters_xy(line_2[0], line_2[1])
    D=par_1[0]*par_2[1]-par_1[1]*par_2[0]
    Dx=par_1[2]*par_2[1]-par_1[1]*par_2[2]
    Dy=par_1[0]*par_2[2]-par_1[2]*par_2[0]
    if D!=0:
        x=Dx/D
        y=Dy/D
        return x, y, 0.0
    else:
        return None

def segments_intersection_xy(line_1, line_2):
    """
    returns the intersection point (coordinates) of two segments 
    """
    int_pt=lines_intersection_xy(line_1, line_2)
    if int_pt!=None:
        len_1=distance_between_two_points(line_1[0], line_1[1])
        len_2=distance_between_two_points(line_2[0], line_2[1])
        dis_11=distance_between_two_points(int_pt, line_1[0])
        dis_12=distance_between_two_points(int_pt, line_1[1])
        dis_21=distance_between_two_points(int_pt, line_2[0])
        dis_22=distance_between_two_points(int_pt, line_2[1])
        if dis_11<=len_1 and dis_12<=len_1 and dis_21<=len_2 and dis_22<=len_2:
            return int_pt
    else:
        return int_pt

def leaf_edges(network):
    """
    returns leaf edges of the "compas" network as a list 
    """ 
    leaf_ver_lis=network.leaves() 
    leaf_edg_lis=[]
    for key in leaf_ver_lis:
        edg=network.vertex_connected_edges(key) 
        leaf_edg_lis.append(edg[0])
    
    return leaf_edg_lis

def leaf_pair_dict(network):
    """
    returns the dict. of leaf pairs with a common vertex (vc): {vc: [(v1,vc), (v2,vc)]
    network: "compas" network
    """
    leaf_edg_lis=leaf_edges(network)
    leaf_ver_lis=network.leaves()
    leaf_nbr_set=set([network.vertex_neighbors(key)[0] for key in leaf_ver_lis])
    leaf_pair_dic={k:[] for k in leaf_nbr_set}
    for k in leaf_pair_dic:
        for edg in leaf_edg_lis:
            if k in edg:
                leaf_pair_dic[k].append(edg)
    
    return leaf_pair_dic
    
def bar_properties(bars, ver_coor):
    """
    returns connectivity matirx, bar lenghts and bar direction cosines (coor_diff/len)
    bars: the paired list of lists of vertices
    ver_coor: vertex coordinates of the bars 
    """
    C=connectivity_matrix(bars, rtype='csr')
    xy=ver_coor.copy()
    uv=C.dot(xy)  # coordinate difference
    bars_len=np.sqrt(np.sum(np.square(uv), axis=1))
    norm_dir=uv/bars_len[:, np.newaxis]  # [[sin_1 cos_1],[sin_2 cos_2]...]
        
    return C, bars_len, np.round_(norm_dir, 2)

def plot_mesh(mesh):
    """
    Compas Mesh Plotter
    """
    all_keys=list(mesh.vertices())
    plotter=MeshPlotter(mesh)
    plotter.draw_vertices(text='key', keys=all_keys, radius=0.01)
    plotter.draw_faces()
    plotter.show()

def plot_network(network):
    """
    Compas Network Plotter
    """
    plotter=NetworkPlotter(network)
    plotter.draw_vertices(radius=0.001, text='key', fontsize=15.0, facecolor=(0, 0, 0)) 
    plotter.draw_edges() 
    plotter.show()

"############################################### AGS ######################################################"

def dual_connectivity_matrix(dic_he, dic_fc, edg_dic):
    """
    returns the connectivity matrix for the dual network
    """
    cs_inc=np.zeros(shape=(len(dic_fc), len(edg_dic)))
    for ver1, dic in dic_he.items():
        for ver2 in dic.keys():
            if dic_he[ver1][ver2] is not None and dic_he[ver2][ver1] is not None:
                if (ver1, ver2) in edg_dic.values():
                    cs_inc[dic_he[ver1][ver2]-1][list(edg_dic.values()).index((ver1, ver2))]=1.0
                else:
                    cs_inc[dic_he[ver1][ver2]-1][list(edg_dic.values()).index((ver2, ver1))]=-1.0
    return cs_inc  

def dual_edges(cs_inc):                   
    """
    returns the connectivity (edge) list of vertices in dual for plotting
    """    
    dual_edg_lis=[]
    for j in range(cs_inc.shape[1]):
        lis_ver=[]
        for i in range(cs_inc.shape[0]):
            if abs(cs_inc[i][j])==1:
                lis_ver.append(i)
        dual_edg_lis.append(lis_ver)
    
    return dual_edg_lis

def dual_coordinates(uv, cs_inc, q_c):
    """
    returns the coordinates of the force (dual) diagram
    """ 
    cs_inc_t=cs_inc.T
    # pick a known point to start force polygon:
    xy_0=np.array([[0.0], [0.0]]).T
    ind_d=0  # index of the known x-y coordinate
    # remove extra row of C* and automatically extra col of C*_tranpose, then find Laplacian
    cs_inc_d=np.delete(cs_inc, ind_d, 0)
    cs_inc_d_t=cs_inc_d.T
    l_s=np.dot(cs_inc_d, cs_inc_d_t)  # Laplacian
    # get the x*, y* dual coordinates
    term_1=np.dot(cs_inc_t[:, [ind_d]], xy_0)
    q_diag=np.diag(q_c.transpose().tolist()[0])  # Q matrix
    term_2=np.dot(q_diag, uv)-term_1
    term_3=np.dot(cs_inc_d, term_2)
    xy_s_d=np.dot(np.linalg.inv(l_s), term_3)
    xy_s=np.vstack((xy_0, xy_s_d))
    
    return xy_s

def rotate_leaves_for_face_rec(ags_net, gtopt_net, plygn, plyln):
    """
    rotates the leaves of the gtopt_net and adds them to ags_net to facilitate face recognition.
    """
    ags_net_rot=gtopt_net.copy()
    leaf_pair_dic=leaf_pair_dict(ags_net)
    for key, pair in leaf_pair_dic.items():
        coor_key=ags_net.vertex_coordinates(pair[0][0])  # the common coordinate
        if len(pair)>1:  # if they are two pairs (sup/load) at one node
            coor_12=ags_net.vertex_coordinates(pair[0][1])
            coor_22=ags_net.vertex_coordinates(pair[1][1])
            plyln_bln_1=is_point_on_polyline(coor_12, plyln.points, tol=0.1)
            plyln_bln_2=is_point_on_polyline(coor_22, plyln.points, tol=0.1)
            if plyln_bln_1 or plyln_bln_2:  # the case when one is on polyline        
                if plyln_bln_1:
                    key_g=pair[0][1]
                    key_o=pair[1][1]
                    coor_g=coor_12
                    coor_o=coor_22
                elif plyln_bln_2:
                    key_g=pair[1][1]
                    key_o=pair[0][1]
                    coor_g=coor_22
                    coor_o=coor_12
                add_vec=add_vectors(vector_create(coor_o, coor_key).tolist(), vector_create(coor_g, coor_key).tolist())
                add_pt=subtract_vectors(coor_key, add_vec)  # bcs the origin is the key_coor
                ags_net_rot.add_vertex(key_g, {'x': add_pt[0], 'y': add_pt[1], 'z': add_pt[2]})
                ags_net_rot.add_vertex(key_o, {'x': coor_o[0], 'y': coor_o[1], 'z': coor_o[2]})
                ags_net_rot.add_edge(key, key_g)
                ags_net_rot.add_edge(key, key_o)
            else:  # the case when both are not on polyline 
                ags_net_rot.add_vertex(pair[0][1], {'x': coor_12[0], 'y': coor_12[1], 'z': coor_12[2]})
                ags_net_rot.add_vertex(pair[1][1], {'x': coor_22[0], 'y': coor_22[1], 'z': coor_22[2]})
                ags_net_rot.add_edge(key, pair[0][1])
                ags_net_rot.add_edge(key, pair[1][1])  
        else:  # for single leaf
            coor_12=ags_net.vertex_coordinates(pair[0][1])
            plyln_bln=is_point_on_polyline(coor_12, plyln.points, tol=0.1)
            if plyln_bln:
                uv=unit_vector(vector_create(coor_key, coor_12))
                if uv[0]-0.0<0.01:  # x=0
                    coor_g=add_vectors(coor_12, (1.0, 0.0, 0.0))
                    plygn_bln=is_point_in_polygon_xy(coor_g, plygn.points)
                    if plygn_bln:
                        coor_g=add_vectors(coor_12, (-1.0, 0.0, 0.0))
                elif uv[1]-0.0<0.01:  # y=0
                    coor_g=add_vectors(coor_12, (0.0, 1.0, 0.0))
                    plygn_bln=is_point_in_polygon_xy(coor_12, plygn.points)
                    if plygn_bln:
                        coor_g=add_vectors(coor_g, (0.0,-1.0, 0.0))            
                ags_net_rot.add_vertex(pair[0][1], {'x': coor_g[0], 'y': coor_g[1], 'z': coor_g[2]})
                ags_net_rot.add_edge(key, pair[0][1])
            else:  # when already in the correct position
                ags_net_rot.add_vertex(pair[0][1], {'x': coor_12[0], 'y': coor_12[1], 'z': coor_12[2]})
                ags_net_rot.add_edge(key, pair[0][1])
    # plot_network(ags_net_rot)
    
    return ags_net_rot

def leaf_edge_dict(edg_dic, network):
    """
    makes a edg_dic of leaf edges
    """
    leaf_ind_edg_dic={}
    leaf_ver_lis=network.leaves()
    for ind, edg in edg_dic.items():
        if edg[0] in leaf_ver_lis or edg[1] in leaf_ver_lis:
            leaf_ind_edg_dic[ind]=edg
  
    return leaf_ind_edg_dic

def get_halfedge_face(network):
    """
    returns halfedge and face dictionary of the network
    """
    net=network.copy()
    mesh=Mesh()
    for key, attr in network.vertices(True):
        mesh.add_vertex(key, x=attr['x'], y=attr['y'], z=attr['z'])
    mesh.halfedge=net.halfedge
    network_find_faces(mesh, mesh.leaves())
    dic_he=mesh.halfedge 
    dic_fc=mesh.face 
    
    return dic_he, dic_fc

def inner_vertices(ver_dic, edg_dic):
    """
    makes a list of inner vertex indeces
    """
    inner=[]
    for ind in ver_dic:
        if len([k for (k, v) in edg_dic.items() if v.count(ind)==1])>1:
            inner.append(ind)
    return inner

def make_network(ver_dic, edg_dic):
    """
    make "compas" network using ver_dic and edg_dic
    """
    net=Network()
    for ind, ver in ver_dic.items():
        net.add_vertex(ind, {'x': ver[0], 'y': ver[1], 'z': ver[2]})
    for edg in edg_dic.values():
        net.add_edge(edg[0], edg[1])
    
    return net
 
"############################################ POST-PROCESSING ############################################"
          
def process_aligned_edges(network):
    """
    Post-Processing Function:
    removes the vertex with degree two where the edges are alined 
    removes corresponding edges
    adds a long edge instead
    returns the cured network and add_rmv_edg_dic={new_edg: (old_edg1, old_edg2)}
    """
    pnt_net=network.copy()
    for key in network.vertex:
        # check which vertex has degree 2
        if pnt_net.vertex_degree(key)==2:
            # find its neighbours
            nbr_lis=pnt_net.vertex_neighbors(key)
            # find the angle between two edges
            vec_1=vector_create(pnt_net.vertex_coordinates(nbr_lis[0]), pnt_net.vertex_coordinates(key))
            vec_2=vector_create(pnt_net.vertex_coordinates(nbr_lis[1]), pnt_net.vertex_coordinates(key))
            ang=angle_between(vec_1, vec_2)
            if round(ang, 2)==0.00 or round(ang, 2)==3.14:
                pnt_net.delete_vertex(key)
                pnt_net.add_edge(nbr_lis[0], nbr_lis[1])
    
    return pnt_net

def process_leaves(network, dic_load, all_bln):
    """
    Post-Processing Function:
    removes the leaf vertices(corresponding edges will be removed automatically)
    returns the cured network
    all_bln: True, removes all leaves, False, just the ones along load
    """
    leaf_ver_lis=network.leaves()
    key_removed={}
    for key in leaf_ver_lis:
        if key in dic_load:
            # there is just one neighbour to the leaf vertex
            key_neighbour=network.vertex_neighbors(key)[0]
            key_removed[key]={key_neighbour: dic_load[key]}
            network.delete_vertex(key)
        if all_bln is True:
            network.delete_vertex(key)
    # to be used in "add_vertex_edge_for_support_load"
    key_removed_dic=key_removed.copy()
    
    return network, key_removed_dic

def check_point_with_vertices(pt_coor, network):
    """
    Process_Crossings Function(used in "add_vertex_at_crossings"): 
    Check if a point is already exists in a list of vertices.
    if True: returns the vertex key
    else: returns None
    Note: precision could be adjusted.
    """
    vertices=network.vertices()
    pt_coor_formatted=geometric_key(pt_coor, precision='3f')  # tolerance=1e-1
    ver_coor_lis=[network.vertex_coordinates(ver) for ver in vertices]
    # map geometric keys (formatted xyz) to network vertices keys
    geo_key_lis=[geometric_key(xyz, precision='3f') for xyz in ver_coor_lis]
    ver_key_lis=[key for key in network.vertices()]
    geo_key_dic=dict(zip(geo_key_lis, ver_key_lis))
    if pt_coor_formatted in geo_key_dic:
        ver_key=geo_key_dic.get(pt_coor_formatted)
        return ver_key
    else:
        return None

def add_vertex_at_crossings(network):
    """
    Process Crossings Function: (used in "sort_vertices_on_edge")
    adds vertex at crossings if it is not already there.
    returns dictionary of added vertices on each edge {crossed_edge_tup:on_edg_ver_key_lis}
    """
    edg_crossings_lis=network_find_crossings(network)  # [(edg1_tup,edg2_tup), (( , ) , ( , )), ...]    
    # make a list of edges which have crossings
    cross_edg_set=set()
    for tup in edg_crossings_lis:
        if (tup[0][1], tup[0][0]) not in cross_edg_set:  # not to add reverse repetative edge
            cross_edg_set.add(tup[0])
        if (tup[1][1], tup[1][0]) not in cross_edg_set:
            cross_edg_set.add(tup[1])
    # make a dict of {crossed_edge_tup:on_edg_ver_key_lis}
    ver_on_edg_dic={edg: [] for edg in cross_edg_set}
    for tup in edg_crossings_lis:
        tup_1=tup[0]
        tup_2=tup[1]
        if tup[0] not in cross_edg_set:  # to reverse the repetative edge
            tup_1=(tup[0][1], tup[0][0])  
        if tup[1] not in cross_edg_set:  # to reverse the repetative edge
            tup_2=(tup[1][1], tup[1][0]) 
        tup_coor_1=(network.vertex_coordinates(tup_1[0]), network.vertex_coordinates(tup_1[1]))
        tup_coor_2=(network.vertex_coordinates(tup_2[0]), network.vertex_coordinates(tup_2[1]))
        x, y, z=segments_intersection_xy(tup_coor_1, tup_coor_2)
        search_key=check_point_with_vertices((x, y, z), network)
        if search_key is None:
            ver_key=network.add_vertex(x=x, y=y, z=z)
        else:
            ver_key=search_key
        if ver_key not in ver_on_edg_dic[tup_1]:   
            ver_on_edg_dic[tup_1].append(ver_key)
        if ver_key not in ver_on_edg_dic[tup_2]: 
            ver_on_edg_dic[tup_2].append(ver_key) 
            
    return ver_on_edg_dic

def sort_vertices_on_edge(network):
    """
    Process_Crossings Function: 
    sorts vertex keys on an edge based on their distance to the edge[0]
    returns dictionary of edge {crossed_edge_tup:sorted_vertiecson_edg_lis} e.g. {(sv,ev):[sv,22,30,51,ev]}
    """
    ver_on_edg_dic=add_vertex_at_crossings(network)
    sorted_ver_on_edg_dic={edg: [] for edg in ver_on_edg_dic}
    for edg, ver_lis in ver_on_edg_dic.items():
        if len(ver_lis)>1:  # case of several crossings on one edge
            dist_dic={}
            dist_lis=[]
            for key in ver_lis:
                edg_0_coor=network.vertex_coordinates(edg[0])
                ver_coor=network.vertex_coordinates(key)
                dist=round(distance_between_two_points(edg_0_coor, ver_coor), 2)
                dist_lis.append(dist)
                dist_dic[str(dist)]=key
            sorted_dist_lis=sorted(dist_lis)
            sorted_ver_lis=[]
            sorted_ver_lis.append(edg[0])
            for dist in sorted_dist_lis:
                sorted_ver_lis.append(dist_dic[str(dist)])
            sorted_ver_lis.append(edg[1])
            sorted_ver_on_edg_dic[edg]=sorted_ver_lis
        else:  # case of one crossing on an edge
            sorted_ver_on_edg_dic[edg]=[edg[0], ver_lis[0], edg[1]]
    
    return sorted_ver_on_edg_dic

def process_crossings(network):
    """
    Post-Processing Function:
    when embedding is unsuccessful, vertices will be added at the crossings.
    the crossed edges will be removed, new edges will be added
    returns the cured network
    """
    # add vertices at crossings and sort the vertices on an edge
    sorted_ver_on_edg_dic=sort_vertices_on_edge(network)
    # delete old crossing edges and add new edges
    for edg, ver_lis in sorted_ver_on_edg_dic.items():
        network.delete_edge(edg[0], edg[1])
        for ind, key_1 in enumerate(ver_lis[:-1]):
            key_2=ver_lis[ind+1]
            network.add_edge(key_1, key_2)
    
    return network

def add_vertex_edge_for_load_support(network, sup_dic, load_dic, bars_len, key_removed_dic):
    """
    Post-Processing Function:
    Adds vertices and edges in accordance with supports and loads
    returns the cured network
    """
    if not key_removed_dic:
        load_sup_dic=merge_two_dicts(sup_dic, load_dic)
    else:
        load_dic_2=load_dic.copy()
        for key in key_removed_dic:
            load_dic_2.pop(key)
            load_dic_2=merge_two_dicts(load_dic_2, key_removed_dic[key])
        load_sup_dic=merge_two_dicts(sup_dic, load_dic_2)
    # define arbitrsry r to be added to get leaf vertex coordinates
    max_len=max(bars_len)
    r=max_len/3.0
    
    # make a polygon and polyline from outer vertices of network
    face_net=Mesh.from_data(network.to_data())
    face_net.halfedge=network.halfedge
    network_find_faces(face_net)
    if 0 in face_net.face and len(face_net.face)>1:
        face_net.delete_face(0)
    if len(face_net.face)==1:  
        ver_lis=[key for key in face_net.vertices()]
    else:
        ver_lis=face_net.vertices_on_boundary(ordered=True)
     
    ver_lis_plyln=ver_lis[:]
    ver_lis_plyln.append(ver_lis[0])
    pt_lis_plygn=[face_net.vertex_coordinates(key) for key in ver_lis]
    pt_lis_plyln=[face_net.vertex_coordinates(key) for key in ver_lis_plyln]
    plygn=Polygon(pt_lis_plygn)
    plyln=Polyline(pt_lis_plyln)

    # add leaf vertices
    for key in load_sup_dic:
        if load_sup_dic[key][0]!=0.0:
            pt_1=add_vectors(network.vertex_coordinates(key), (+r, 0.0, 0.0))
            plyln_bln=is_point_on_polyline(pt_1, plyln.points, tol=0.001)
            plygn_bln=is_point_in_polygon_xy(pt_1, plygn.points)
            if plyln_bln or plygn_bln:
                pt_1=add_vectors(network.vertex_coordinates(key), (-r, 0.0, 0.0))
            key_2=network.add_vertex(x=np.asscalar(pt_1[0]), y=pt_1[1], z=0.0)
            network.add_edge(key, key_2)
        if load_sup_dic[key][1]!=0.0:
            pt_2=add_vectors(network.vertex_coordinates(key), (0.0,+r, 0.0))
            plyln_bln=is_point_on_polyline(pt_2, plyln.points, tol=0.001)
            plygn_bln=is_point_in_polygon_xy(pt_2, plygn.points)
            if plyln_bln or plygn_bln:
                pt_2=add_vectors(network.vertex_coordinates(key), (0.0,-r, 0.0))
            key_2=network.add_vertex(x=pt_2[0], y=np.asscalar(pt_2[1]), z=0.0)
            network.add_edge(key, key_2)
    
    return network, plygn, plyln

def ags_inputs(network):
    """
    Post-Processing Function:
    saves edg_dic_GT and ver_dic_GT files to be used in AGS
    """
    ags_net=Network()  # make a new network 
    new_key_lis=range(len(network.vertex))
    old_key_lis=list(network.vertices())
    map_key_dic=dict(zip(old_key_lis, new_key_lis))  # {old_key:new_key}
 
    for ver in map_key_dic:
        ver_coor=network.vertex
        ags_net.add_vertex(key=map_key_dic[ver], x=ver_coor[ver]['x'], y=ver_coor[ver]['y'], z=ver_coor[ver]['z'])
    for edg in network.edges():
        u=map_key_dic[edg[0]]
        v=map_key_dic[edg[1]]
        if (u, v) not in ags_net.edges() and (v, u) not in ags_net.edges():
            ags_net.add_edge(u, v)
    
    ver_dic={}
    edg_dic={}
    for key in ags_net.vertex:
        ver_dic[key]=ags_net.vertex_coordinates(key)
    for ind_edg, edg in enumerate(ags_net.edges()):
        edg_dic[ind_edg]=edg  
    
    with open('ver_dic_GT.p', 'wb') as fp: 
        pickle.dump(ver_dic, fp, protocol=2)
    with open('edg_dic_GT.p', 'wb') as fp:
        pickle.dump(edg_dic, fp, protocol=2)
        
    return ags_net        

"######################################## PROCESS-DUAL-DIAGRAMS ###########################################"

def map_vertices_with_similar_coordinates(force_orig_net):
    """
    returns a mapping dictionary of vertices with similar coordinates
    (i.e. maps all similar coordinates' keys to a single key
    e.g. when key0-2 have the same coordinates: {key0:key0, key1:key0, key2:key0})
    used in "update_dual_mapping_1"
    """
    # make a dict of {coor: [list of keys with similar coors]}
    map_coor_key_dic={}
    for key in force_orig_net.vertex:
        coor=force_orig_net.vertex_coordinates(key)
        coor_formatted=geometric_key(coor, precision='0f')  # abut_rec "d"
        if coor_formatted not in map_coor_key_dic:
            map_coor_key_dic[coor_formatted]=[]
        map_coor_key_dic[coor_formatted].append(key)
    # map list elements of map_coor_key_dic to the first element of each list
    map_key_dic={}
    for lis in map_coor_key_dic.values():
        for key in lis:
            map_key_dic[key]=lis[0]

    return map_key_dic    

def update_dual_mapping_1(force_orig_net, map_edg_orig_dic):
    """
    updates map_edg_dic_orig by removing the repetative vertices in dual
    returns the intermediate temporary mapping dict
    """
    map_key_dic=map_vertices_with_similar_coordinates(force_orig_net)
    # update duals edg mapping
    map_edg_temp_dic={}
    for edg, dual_edg in map_edg_orig_dic.items():
        new_dual_edg=[map_key_dic[dual_edg[0]], map_key_dic[dual_edg[1]]]
        if [new_dual_edg[1], new_dual_edg[0]] in map_edg_temp_dic.values():
            map_edg_temp_dic[edg]=(new_dual_edg[1], new_dual_edg[0])
        else:
            map_edg_temp_dic[edg]=tuple(new_dual_edg)

    return map_edg_temp_dic

def outer_inner_faces(fc_dic, network):
    """
    returns the dictionaries of outer and inner faces
    used in "update_dual_mapping_2"
    """
    leaf_set=set(network.leaves())
   
    out_fc_set=set()
    for fc_key, key_lis in fc_dic.items():
        for leaf in leaf_set:
            if leaf in key_lis:
                out_fc_set.add(fc_key)
    
    out_fc_dic={key: fc_dic[key] for key in out_fc_set}
    in_fc_dic={key: fc_dic[key] for key in fc_dic.keys() if key not in out_fc_set}

    return out_fc_dic, in_fc_dic

def make_network_from_face(network, key_lis):
    """
    makes a network from list of vertex keys of a face
    used in "find_inner_face_corners" and "find_outer_face_corners"
    """
    net=Network()
    for key in key_lis:
        xyz=(network.vertex[key]['x'], network.vertex[key]['y'], network.vertex[key]['z'])
        net.add_vertex(key, {'x': xyz[0], 'y': xyz[1], 'z': xyz[2]})   
    for edg in network.edges():
        if edg[0] in key_lis and edg[1] in key_lis:
            net.add_edge(edg[0], edg[1])     
    
    return net
    
def find_inner_face_corners(network, in_dic_fc): 
    """
    finds the corners of an inner face (i.e. degree 2 w/ non-180 deg angles), used in "update_dual_mapping_2"
    in_fc_dic: inner face dictionary from the FaceNetwork of the network
    returns: corner_dic = {face_key : list of corner vertices}
    """
    corner_dic={k: [] for k in in_dic_fc}
    in_corner_set=set()
    for fc_key, lis in in_dic_fc.items():    
        # make a new newtork for each face
        net=make_network_from_face(network, lis)
        # now find the corners by looking at the angles
        for key in lis:
            nbr_lis=net.vertex_neighbors(key)
            if len(nbr_lis)>1:
                vec_1=vector_create(net.vertex_coordinates(nbr_lis[0]), net.vertex_coordinates(key))
                vec_2=vector_create(net.vertex_coordinates(nbr_lis[1]), net.vertex_coordinates(key))
                ang=angle_between(vec_1, vec_2)
                if round(ang, 2)!=0.00 and round(ang, 2)!=3.14:
                    corner_dic[fc_key].append(key)
                    in_corner_set.add(key)
                elif network.vertex_degree(key)>2:
                    corner_dic[fc_key].append(key)
                    in_corner_set.add(key)    

    return corner_dic, in_corner_set

def find_outer_face_corners(network, out_dic_fc, in_corner_set): 
    """
    finds the corners of an outer face(i.e degree 2 w/ non-180 deg angles), used in "update_dual_mapping_2"
    out_fc_dic: outer face dictionary from the FaceNetwork of the network
    in_corner_lis: list of corners for the inner faces
    returns: corner_dic = {face_key : list of corner vertices}
    """
    corner_dic={k: [] for k in out_dic_fc}
    for fc_key, lis in out_dic_fc.items():
        # make a new newtork for each face
        net=make_network_from_face(network, lis)
        # now find the corners by looking at the angles
        for key in lis:
            nbr_lis=net.vertex_neighbors(key)
            if len(nbr_lis)==2:  # non-leaf vertices
                vec_1=vector_create(net.vertex_coordinates(nbr_lis[0]), net.vertex_coordinates(key))
                vec_2=vector_create(net.vertex_coordinates(nbr_lis[1]), net.vertex_coordinates(key))
                ang=angle_between(vec_1, vec_2)
                if round(ang, 2)!=0.00 and round(ang, 2)!=3.14:
                    corner_dic[fc_key].append(key)
                elif key in in_corner_set:  # if it is a corner of inner face, should be the corner of outer face
                    corner_dic[fc_key].append(key)
            elif len(nbr_lis)==1:  # leaf vertices
                corner_dic[fc_key].append(key)        
 
    return corner_dic 

def find_face_aligned_edges(ver_lis, corner_lis):
    """
    finds the aligned edges of a face, used in "update_dual_mapping_2"
    ver_lis= vertice list of a face
    corner_lis= vertex corners of a face
    returns: face_edg_dic = {face_corner_edg : list of vertices on the long edge}
    """
    ind_mod=len(ver_lis)
    face_edg_dic={}
    # start from a corner and go forward to get to the next corner 
    for key in corner_lis:
            align_ver_lis=[]
            align_ver_lis.append(key)  # append the 1st corner
            cor_ind=ver_lis.index(key)
            ind=1
            while ver_lis[(cor_ind+ind)%ind_mod] not in corner_lis:
                align_ver_lis.append(ver_lis[(cor_ind+ind)%ind_mod])
                ind+=1
            align_ver_lis.append(ver_lis[(cor_ind+ind)%ind_mod])  # append the 2nd corner         
            face_edg_dic[(align_ver_lis[0], align_ver_lis[-1])]=align_ver_lis
                  
    return face_edg_dic            

def update_dual_mapping_2(form_orig_net, map_edg_temp_dic, old_edg_f_dic):
    """
    updates map_edg_temp_dic, maps aligned short edges to one long edge
    returns the final mapping dict
    """
    _, fc_dic=get_halfedge_face(form_orig_net)
    out_fc_dic, in_fc_dic=outer_inner_faces(fc_dic, form_orig_net)
    
    # get the dictionary of corners according to each face_key
    in_corner_dic, in_corner_set=find_inner_face_corners(form_orig_net, in_fc_dic)
    out_corner_dic=find_outer_face_corners(form_orig_net, out_fc_dic, in_corner_set)
    corner_dic=merge_two_dicts(in_corner_dic, out_corner_dic)

    # get aligned vertices from one corner to the next corner (each face at a time)
    # map the corner to corner (new added edges) to one of the duals
    map_edg_dic={}
    new_edg_f_dic={} 
    for fc_key, lis in fc_dic.items():
        face_edg_dic=find_face_aligned_edges(lis, corner_dic[fc_key])  # {(new_edg):[new_edg_sp, old_key_1, old_key_2, ..., new_edg_ep]}
        for new_edg in face_edg_dic:
            if new_edg not in map_edg_dic and (new_edg[1], new_edg[0]) not in map_edg_dic:
                rmv_edg=(face_edg_dic[new_edg][0], face_edg_dic[new_edg][1])  # to get the dual edge from temp map (new_edg_sp, old_key_1)
                if rmv_edg in map_edg_temp_dic:
                    map_edg_dic[new_edg]=tuple(map_edg_temp_dic[rmv_edg])           
                    # map forces to the new edges
                    new_edg_f_dic[new_edg]=old_edg_f_dic[rmv_edg]
                elif (rmv_edg[1], rmv_edg[0]) in map_edg_temp_dic:
                    rmv_edg=(rmv_edg[1], rmv_edg[0])                 
                    map_edg_dic[new_edg]=tuple(map_edg_temp_dic[rmv_edg])
                    # map forces to the new edges
                    new_edg_f_dic[new_edg]=old_edg_f_dic[rmv_edg]
    
    return map_edg_dic, new_edg_f_dic

def make_new_network (orig_net, edg_lis):
    """
    makes a new newtork according to new edges
    """
    new_net=Network()
    for edg in edg_lis:
        coor0=orig_net.vertex_coordinates(edg[0])
        coor1=orig_net.vertex_coordinates(edg[1])
        if edg[0] not in new_net.vertex:
            new_net.add_vertex(edg[0], {'x': coor0[0], 'y': coor0[1], 'z': coor0[2]})
        if edg[1] not in new_net.vertex:    
            new_net.add_vertex(edg[1], {'x': coor1[0], 'y': coor1[1], 'z': coor1[2]})
        new_net.add_edge(edg[0], edg[1])

    return new_net
        
def rotate_dual(force_net, ANG):
    """
    rotates the force_net coordinates 90 deg counterclockwise
    """
    # save the rotated verteces in ver_coor_90_dic    
    ver_coor_90_dic={}
    AX=[0.0, 0.0, 1.0]  # normal to the plane of rotation
    ORG=[0.0, 0.0, 0.0]
    for key in force_net.vertex:
        coor=force_net.vertex_coordinates(key)
        pt=rotate_points([coor], ANG, AX, ORG) 
        ver_coor_90_dic[key]=np.round_(pt[0], 3).tolist()            
    # make new rotated dual network 
    force_90_net=Network()
    for key, coor in ver_coor_90_dic.items():
        force_90_net.add_vertex(key, {'x': coor[0], 'y': coor[1], 'z': coor[2]})       
    for edg in force_net.edges():
        force_90_net.add_edge(edg[0], edg[1])    

    return force_90_net

"############################# STRESS FIELDS ########################################"

def sf_cl_anchor_lines(sc_pts_lis, pts_lis, CONST):
    """
    finds the anchorage lines and stress field center lines in sf (plotting purposes)
    returns the list of lines, 3 lines for each case: [(pt1_CL),(pt2_CL)],[(pt1_AL1),(pt2_AL1)] ,[(pt1_AL2),(pt2_AL2)] ])
    sc_pts_lis: scaled polygon (dual) points (list of two points)
    pts_lis: truss member points (list of two points)
    CONST: 0, 1 or 2 refrains from producing anchors if the vertex is leaf (see "constant_stress_fields")
    """
    mid_pt=middle_point(sc_pts_lis[0], sc_pts_lis[1])
    line_pts_lis=[]
    line_pts=[]
    # for ties 
    for pt in pts_lis:
        new_mid_pt=add_vectors(mid_pt, pt)
        line_pts.append((new_mid_pt[0], new_mid_pt[1]))
    line_pts_lis.append(line_pts)
    # for anchor plates 
    if CONST==0:
        pts_lis=[pts_lis[1]]
    if CONST==1:
        pts_lis=[pts_lis[0]]
    for pt in pts_lis:  
        line_pts=[]
        for dl_pt in sc_pts_lis:            
            new_pt=add_vectors(dl_pt, pt)
            line_pts.append((new_pt[0], new_pt[1]))
        line_pts_lis.append(line_pts)
    
    return line_pts_lis  # consists of point pairs (2 or 3 lines: 1 tie, 1 or 2 anchor lines)

def minkowski_sum(sc_pts_lis, pts_lis):
    """
    Just compression cases!
    calculates a list of sum of all the points in the two input lists (Minkowski Sum)
    sc_pts_lis: scaled polygon (dual) points (list of two points)
    pts_lis: structure element points (list of two points) 
    returns the scipy convex hull of the new points
    """
    hull_pts_lis=[]
    for du_pt in sc_pts_lis:
        for pt in pts_lis:
            new_pt=add_vectors(du_pt, pt)
            hull_pts_lis.append([new_pt[0], new_pt[1]])
    pts=np.array(hull_pts_lis)  # consist of 4 poitns
    hull=ConvexHull(pts)  # scipy   
    
    return hull

def plot_stress_fields(hull_lis, ten_lines_dic):
    """
    !!! "unadjusted" stress fields !!!
    plots the convex_hulls of the summed polylines/polygons for compression case
    plots the element CL and anchorage areas for tension cases
    """
    fig=plt.figure(num='Stress Fields')
    ax=fig.gca()
    ax.set_title('stress fields', fontsize=15)
    # Compression struts
    for hull in hull_lis:
        line_segments=[hull.points[simplex] for simplex in hull.simplices]
        ax.add_collection(LineCollection(line_segments, colors='k', linestyle='solid', linewidths=1.0, zorder=1))
    # Tention ties and anchoring areas
    for lines in ten_lines_dic.values():
        x, y=[lines[0][0][0], lines[0][1][0]] , [lines[0][0][1], lines[0][1][1]]
        ax.plot(x, y, color='r', zorder=4)  # tension tie
        ax.add_collection(LineCollection(lines[1:], colors='k', linestyle='solid', linewidths=3.5, zorder=4))  # anchor plates
    ax.axis('off')
    ax.plot()
    
def sf_rhino_inputs(hull_lis, ten_lines_dic):
    """
    create hull_dic to save as rhino input
    hull_dic = {hull_ind = [seg1[[pt1][pt2]], seg2[[][]], ...], ...}
    saves lines_dic
    """
    hull_dic={ind: [] for ind in range(len(hull_lis))}
    for ind, hull in enumerate(hull_lis):
        for simplex in hull.simplices:
            hull_dic[ind].append(hull.points[simplex].tolist())     

    with open('hull_dic.p', 'wb') as fp:
        pickle.dump(hull_dic, fp, protocol=2)
    with open('ten_lines_dic.p', 'wb') as fp:
        pickle.dump(ten_lines_dic, fp, protocol=2)  
