'''
Created on 06.07.2018, updated on 24.02.2020

@author: Salma Mozaffari, ETH Zurich

Comment: Integration of AGS and LAYOPT:
        - Calculates and plots form and force diagrams
        - Calculates and plots the constant stress fields

Before running the code:
        - install required packages: matplotlib, numpy, scipy, cvxpy, and compas (has tested with compas version 15.4).
        - in the RUN section, update the directory to the "MeshObjects" folder and the specific mesh object according to your directory.
        - in the INPUTS section, uncomment the "dic_sup" and "dic_load" to impose boundary conditions (supports and external loads) on the mesh object.
        - in the RUN section, uncomment "hf.plot_network(AO.dic_attr['gt_net'])" to see/check ground truss and the nodal locations of supports/external loads

After a correct run:
        - three plots of form diagram (optimized truss), force diagram, and uneditted stress fields will appear.
        - the necessary ".p" files will be saved in "Source" folder to be used for drawing in Rhino using "Draw_Diagrams_Rh" and "Draw_SF_Rh".
'''
"################################################ IMPORTS ###################################################"

import Helper_Funcs_2D as hf
import matplotlib.pyplot as plt
import itertools
import cvxpy as cvx
import numpy as np
import pickle
from scipy.sparse import csc_matrix
from compas.datastructures import Mesh
from compas.datastructures.network import Network
from compas.numerical.matrices import connectivity_matrix
from compas.numerical.matrices import equilibrium_matrix
from compas.geometry import scale_points
import os
BASEDIR = os.path.dirname(os.path.realpath(__file__))

"################################################ CLASS ###################################################"

class AGS_OPT_2D (object):

    def __init__(self, dic_sup, dic_load):
        """
        dic_sup: dictionary of supports
                {mesh_node_#: [fix direction in x, fix direction in y]}, list elements: 0 or 1
        dic_load: dictionary of external loads
                {mesh_node_#: [external load in x direction, load in y direction]}
        """
        self.sig_t=4350.0  # steel factored yield stress: 435 MPa [N/mm^2], 435000 [kN/m^2], 4350 [(kN/100)/m^2]
        self.sig_c=130.0  # concrete factored yield stress: 13 MPa [N/mm^2]
        self.den_tol=0.1  # tolerance on force density
        self.load=dic_load  #loads are in [kN/100]
        self.sup=dic_sup
        self.dic_attr={}  # initiate the dictionary of class sttributes

    def create_mesh(self, filepath):
        """
        creates a "compas" mesh from a Rhino mesh object (saved as .obj file)
        filepath: file path in string format
        saves: compas mesh object
        """
        mesh=Mesh.from_obj(filepath)

        self.dic_attr['mesh']=mesh

    def struct_domain(self, mesh):
        """
        mesh: Mesh datastrcuture
        """
        ver_coor=np.array([mesh.vertex_coordinates(key, 'xy') for key in mesh.vertex])
        elem_dic=mesh.face

        self.dic_attr['ver_coor']=ver_coor
        self.dic_attr['elem_dic']=elem_dic

    def generate_ground_truss(self):
        """
        generate ground truss element properties
        """
        ver_coor=self.dic_attr['ver_coor']
        elem_dic=self.dic_attr['elem_dic']

        bars=set()
        for elm in list(list(elem_dic.values())):
            combis=set(itertools.combinations(sorted(elm), 2))
            bars=bars.union(combis)

        C, bars_len, norm_dir=hf.bar_properties(bars, ver_coor)

        self.dic_attr['bars']=bars
        self.dic_attr['bars_len']=bars_len
        self.dic_attr['norm_dir']=norm_dir
        self.dic_attr['conn_mat']=C

    def network_from_ground_truss(self):
        """
        generate compas network of ground truss and saves it in dic_attr['gt_net']
        (the vertices are the same as the mesh vertices)
        """
        mesh=self.dic_attr['mesh']
        bars=self.dic_attr['bars']
        gt_net=Network()
        for v_key, v_attr in mesh.vertex.items():
            gt_net.add_node(v_key, v_attr)
        for edg in bars:
            gt_net.add_edge(edg[0], edg[1])

        self.dic_attr['gt_net']=gt_net

    def post_processing(self, edg_bln, leaf_bln, cross_bln, add_bln, ags_bln, network):
        """
        post-processes the ground/optimized truss to be used in algebraic graph statics (making planar, adding leaves, ...)
        edg_bln= True, when removing aligned edges are desired.
        leaf_bln= True, when removing leaf edges and vertices is desired.
        cross_bln = True, when add a vertex at the crossings elements of the ground truss
        add_bln= True, when adding edges and vertices for sup/load is desired
        ags_bln: True, when need to save edg_dic and ver_dic to be used in AGS
        network: the network that is going to be post-procesed.
        !saves ver_dic_GT and edg_dic_GT files, inside 'ags_inputs'!
        """
        ppnt_net=network.copy()

        # remove alined edges, replace with a longer edge
        if edg_bln is True:
            ppnt_net=hf.process_aligned_edges(ppnt_net)

        # remove leaf vertices and their keys
        if leaf_bln is True:
            ppnt_net, key_removed_dic=hf.process_leaves(ppnt_net, self.load, True)
        else:
            key_removed_dic={}

        # add vertices at crossings
        if cross_bln:
            ppnt_net=hf.process_crossings(ppnt_net)

        # add edges and vertices for supports and external loads
        if add_bln is True:
            ppnt_net, plygn, plyln=hf.add_vertex_edge_for_load_support(ppnt_net, self.sup, self.load, self.dic_attr['bars_len'], key_removed_dic)
            self.dic_attr['plygn']=plygn
            self.dic_attr['plyln']=plyln

        # # save edg_dic and ver_dic to be used in AGS
        if ags_bln is True:
            ppnt_net=hf.ags_inputs(ppnt_net)

        return ppnt_net

    def inc_mat(self):
        """
        constructs incident C matrix for form diagram
        """
        edg_dic=self.dic_attr['edg_dic']
        C=connectivity_matrix(list(edg_dic.values()), rtype='csc')
        c_inc=C.transpose()

        self.dic_attr['c_inc']=c_inc

    def inc_s_mat(self):
        """
        constructs incident C* matrix for force diagram
        """
        # Get face info first
        face_net=hf.rotate_leaves_for_face_rec(self.dic_attr['ags_net'], self.dic_attr['gtopt_net'], self.dic_attr['plygn'], self.dic_attr['plyln'])
        self.dic_attr['dic_he'], self.dic_attr['dic_fc']=hf.get_halfedge_face(face_net)
        cs_inc=hf.dual_connectivity_matrix(self.dic_attr['dic_he'], self.dic_attr['dic_fc'], self.dic_attr['edg_dic'])
        self.dic_attr['cs_inc']=cs_inc

    def coor_dif(self):
        """
        constructs coordinate difference vectors
        """
        ver_dic=self.dic_attr['ver_dic']
        c_inc=self.dic_attr['c_inc']
        xyz_coor=np.array(list(ver_dic.values()))
        uv_coor_diff=c_inc.T*xyz_coor[:, 0:2]  # xyz_coor[:, 0:2] to remove z_coor
        self.dic_attr['uv_dif']=uv_coor_diff

    def equil_mat(self):
        """
        constructs equilibrium matrix A
        """
        ver_dic=self.dic_attr['ver_dic']
        edg_dic=self.dic_attr['edg_dic']
        C=self.dic_attr['c_inc'].T
        xyz=list(ver_dic.values())
        inner=hf.inner_vertices(ver_dic, edg_dic)
        a_equ=equilibrium_matrix(C, xyz, inner, rtype='csc')

        self.dic_attr['inner_ver_lis']=inner
        self.dic_attr['equil_mat']=a_equ

    def fixed_dof(self):
        """
        constructs a set for the fixed directions (i.e. supports)
        """
        n_ver=len(self.dic_attr['inner_ver_lis'])
        sup=self.sup
        fix_dof=set()
        for key in sup:
            if sup[key][0]!=0:
                fix_dof.add(key)
            if sup[key][1]!=0:
                fix_dof.add(key+n_ver)

        self.dic_attr['fix_dof']=fix_dof

    def free_dof(self):
        """
        constructs a list for the non-fixed directions
        """
        fix_dof=self.dic_attr['fix_dof']
        inner=self.dic_attr['inner_ver_lis']
        free_dof=[ind for ind in range(2*len(inner)) if ind not in fix_dof]

        self.dic_attr['free_dof']=free_dof

    def load_vec(self):
        """
        constructs the sparse load vector
        """
        load=self.load
        free_dof=self.dic_attr['free_dof']
        n_ver=len(self.dic_attr['inner_ver_lis'])
        n_load=sum([np.count_nonzero(np.array(load[key])) for key in load])
        data=[]
        rows=[]
        cols=[0]*n_load
        for key, lis in load.items():
            if lis[0]!=0.0:
                data.append(lis[0])
                rows.append(key)
            if lis[1]!=0.0:
                data.append(lis[1])
                rows.append(key+n_ver)
        # complete load_den_vec (ZERO at support dofs)
        f_vec_c=csc_matrix((data, (rows, cols)), shape=(2*n_ver, 1))
        # keep the non-fixed dofs elements in f_vec_c
        f_vec=f_vec_c[list(free_dof), :]

        self.dic_attr['load_dof']=rows
        self.dic_attr['f_vec']=f_vec

    def lin_prog_inputs(self):
        """
        prepares the required input for lin_prog
        """
        edg_dic=self.dic_attr['edg_dic']
        ags_net=self.dic_attr['ags_net']
        equil_mat=self.dic_attr['equil_mat']
        free_dof=self.dic_attr['free_dof']

        leaf_ind_edg_dic=hf.leaf_edge_dict(edg_dic, ags_net)
        # pick one leaf edge and get the length
        leaf_edg=list(leaf_ind_edg_dic.values())[0]
        leaf_len=hf.distance_between_two_points(ags_net.node_coordinates(leaf_edg[0]), ags_net.node_coordinates(leaf_edg[1]))
        # find non-leaf edge indices and save them
        non_leaf_edg_ind_lis=[ind for ind in edg_dic if ind not in list(leaf_ind_edg_dic.keys())]
        # pick the needed bars_len for opt
        ver_coor=np.array([ags_net.node_coordinates(key, 'xy') for key in ags_net.node])
        all_bars=np.array(list(edg_dic.values()))
        _, all_bars_len, all_norm_dir=hf.bar_properties(all_bars, ver_coor)
        bars_len=np.array([edg_len for ind, edg_len in enumerate(all_bars_len) if ind in non_leaf_edg_ind_lis])
        # find straight_diagonal bars for additional opt constraints
        bars=[edg for ind, edg in enumerate(all_bars) if ind in non_leaf_edg_ind_lis]
        norm_dir=np.array([edg_norm for ind, edg_norm in enumerate(all_norm_dir) if ind in non_leaf_edg_ind_lis])
        straight_bars_dic={}
        diag_bars_dic={}

        for ind, edg in enumerate(bars):
            ang=round(np.arctan2(*norm_dir[ind])*180/np.pi, 1)  # angle with y axis
            if ang==0.0 or ang==180.0 or ang==-180.0 or ang==90.0 or ang==-90.0:
                straight_bars_dic[ind]=edg
            else:
                diag_bars_dic[ind]=edg

        # remove the unnecessary rows/cols of the equil_mat
        row_ind=np.array(list(free_dof))
        col_ind=np.array(non_leaf_edg_ind_lis)
        equil_mat_c=equil_mat[:, col_ind]  # Ad of AGS
        equil_mat_cc=equil_mat[row_ind, :][:, col_ind]  # Ad of AGS with the fixed dofs removed

        self.dic_attr['equil_mat_c']=equil_mat_c
        self.dic_attr['equil_mat_cc']=equil_mat_cc
        self.dic_attr['non_leaf_edg_ind_lis']=non_leaf_edg_ind_lis
        self.dic_attr['leaf_len']=leaf_len
        self.dic_attr['bars_len']=bars_len
        self.dic_attr['straight_bars_dic']=straight_bars_dic
        self.dic_attr['diag_bars_dic']=diag_bars_dic

    def lin_prog(self):
        """
        runs the optimization (linaer programming)
        """
        n_bars=len(self.dic_attr['bars_len'])
        bars_len=np.reshape(self.dic_attr['bars_len'], (n_bars,1))
        A_eq=self.dic_attr['equil_mat_cc']
        b_eq=self.dic_attr['f_vec']  # this is the real load (not the density)
        a_var=cvx.Variable((n_bars,1))
        q_var=cvx.Variable((n_bars,1))
        q_var_len=cvx.multiply(bars_len, q_var)  #q_var is the force density
        vol=bars_len.T@a_var
        ind_diag=list(self.dic_attr['diag_bars_dic'].keys())
        constraints=[A_eq@q_var==b_eq, a_var>=10e-8, q_var_len-self.sig_t*a_var<=0.0, -q_var_len-self.sig_c*a_var<=0.0, q_var[ind_diag]<=0.0]
        prob=cvx.Problem(cvx.Minimize(vol), constraints)
        prob.solve()
        print ("status:", prob.status)
        print ("optimized truss volume: [m^3]", round(vol.value[0][0],2))

        self.dic_attr['q_bars']=q_var.value

    def leaf_edg_dof_mapping(self):
        """
        constructs dictionaries:
        leaf edge index and its related degree of freedom: {ind: dof}
        leaf edge index and its related unit vector (from stracture outward direction): {ind: unit vec array}
        """
        edg_dic=self.dic_attr['edg_dic']
        ags_net=self.dic_attr['ags_net']
        fix_dof=list(self.dic_attr['fix_dof'])
        load_dof=self.dic_attr['load_dof']
        n_ver=len(self.dic_attr['inner_ver_lis'])

        leaf_dof=fix_dof+load_dof
        leaf_edg_ind_lis=list(hf.leaf_edge_dict(edg_dic, ags_net).keys())
        leaf_ver_lis=ags_net.leaves()
        leaf_edg_dof_dic={}  # {leaf_edg_ind: dof}
        edg_unit_vec_dic={}  # {leaf_edg_ind: unit_vecotr}  #to be used in "make_q_complete"

        # construct edg_unit_vec_dic
        for ind in leaf_edg_ind_lis:
            edg=edg_dic[ind]
            # find non-leaf vertex (=key or sp) of the leaf edge
            if edg[0] not in leaf_ver_lis:
                sp=edg[0]  # key
                ep=edg[1]  # to get correct unit vec direction in "edg_unit_vec"
            if edg[1] not in leaf_ver_lis:
                sp=edg[1]  # key
                ep=edg[0]  # to get correct unit vec direction in "edg_unit_vec"
            edg_unit_vec=hf.unit_vector(hf.vector_create(ags_net.node_coordinates(sp), ags_net.node_coordinates(ep)))
            edg_ang=round(hf.angle_between(edg_unit_vec, (1.0, 0.0, 0.0)), 2)
            edg_unit_vec_dic[ind]=edg_unit_vec
            # find which dofs has the vertex key (sp) and which leaf edge has which dof
            x_dof=sp
            y_dof=sp+n_ver
            if x_dof in leaf_dof and edg_ang==0.0 or edg_ang==3.14:
                leaf_edg_dof_dic[ind]=x_dof
            elif y_dof in leaf_dof:
                leaf_edg_dof_dic[ind]=y_dof

        self.dic_attr['leaf_edg_dof_dic']=leaf_edg_dof_dic
        self.dic_attr['edg_unit_vec_dic']=edg_unit_vec_dic

    def make_q_complete(self):
        """
        returns the complete force density vector
        """
        edg_dic=self.dic_attr['edg_dic']
        equil_mat_c=self.dic_attr['equil_mat_c']  # complete equil. matrix, still without columns related to leaf edges
        q_bars_arr=self.dic_attr['q_bars']
        leaf_edg_dof_dic=self.dic_attr['leaf_edg_dof_dic']
        non_leaf_edg_ind_lis=self.dic_attr['non_leaf_edg_ind_lis']
        edg_unit_vec_dic=self.dic_attr['edg_unit_vec_dic']
        leaf_len=self.dic_attr['leaf_len']

        f_vec_c=(1.0/leaf_len)*equil_mat_c*q_bars_arr  # complete leaf densities vector (NON-ZERO at support and load dofs)

        q_c=np.zeros((len(edg_dic), 1))
        for ind in edg_dic:
            if ind in leaf_edg_dof_dic:
                non_zero_ind=np.nonzero(edg_unit_vec_dic[ind])  # gets non-zero element of unit vec array
                q_c[ind][0]=f_vec_c[leaf_edg_dof_dic[ind]]*edg_unit_vec_dic[ind][non_zero_ind]
            else:
                q_c[ind][0]=q_bars_arr[non_leaf_edg_ind_lis.index(ind)]

        self.dic_attr['q_c']=q_c

    def dual_coor(self):
        """
        find the the locations of the vertices for dual
        """
        uv_dif=self.dic_attr['uv_dif']
        cs_inc=self.dic_attr['cs_inc']
        q_c=self.dic_attr['q_c']

        xy_s=hf.dual_coordinates(uv_dif, cs_inc, q_c)

        self.dic_attr['xy_s']=xy_s

    def plot_dual_diagrams(self):
        """
        plots form and force diagrams
        """
        ver_dic=self.dic_attr['ver_dic']
        edg_dic=self.dic_attr['edg_dic']
        cs_inc=self.dic_attr['cs_inc']
        xy_s=self.dic_attr['xy_s']
        q=self.dic_attr['q_c']
        # make a connectivity list of vertices in dual for plotting
        dual_edg_lis=hf.dual_edges(cs_inc)

        ### plot force diagram ###
        fig_1=plt.figure(num='Force Diagram')
        ax_1=fig_1.gca()
        ax_1.set_title('force diagram', fontsize=15)
        ax_1.axis('off')
        for ind, edg in enumerate(dual_edg_lis):
            xd, yd=[xy_s[edg[0]][0], xy_s[edg[1]][0]], [xy_s[edg[0]][1], xy_s[edg[1]][1]]
            if q[ind][0]>self.den_tol:  # tension
                ax_1.plot(xd, yd, color='r')
            elif q[ind][0]<-self.den_tol:  # compression
                ax_1.plot(xd, yd, color='b')
           
        ### plot form diagram ###
        fig_2=plt.figure(num='Form Diagram')
        ax_2=fig_2.gca()
        ax_2.set_title('form diagram', fontsize=15)
        ax_2.axis('off')
        for ind, edg in enumerate(list(edg_dic.values())):
            x, y=[ver_dic[edg[0]][0], ver_dic[edg[1]][0]], [ver_dic[edg[0]][1], ver_dic[edg[1]][1]]
            if q[ind][0]>self.den_tol:  # tension
                ax_2.plot(x, y, linewidth=1.0, color='r') 
            elif q[ind][0]<-self.den_tol:  # compression
                ax_2.plot(x, y, linewidth=1.0, color='b')  
            else:
                ax_2.plot(x, y, linewidth=0.1, color='green')
        
        self.dic_attr['dual_edg_lis']=dual_edg_lis

    def update_duals_info(self):
        """
        makes new dictionaries for dual graphs' vertices and edges (i.e. removing zero elements)
        forms a mapping dictinary between dual edges {( , ): [ , ]}
        """
        q_lis=self.dic_attr['q_c'].tolist() # based on dic_attr['edg_dic'] indeces(indeces of original ags_net)
        edg_dic=self.dic_attr['edg_dic']
        ver_dic=self.dic_attr['ver_dic']
        dual_edg_lis=self.dic_attr['dual_edg_lis']  # dual edges
        xy_s=self.dic_attr['xy_s']  # array-dual coordinates
        nr, _=xy_s.shape
        xyz_s=np.c_[xy_s, np.zeros(nr)]
        xyz_s_lis=xyz_s.tolist()

        edg_dic_new={}  # {index: (non-zero_edge)}
        ver_dic_new={}  # {key: [ver_coor]}
        dic_dual_edg={}  # {index: [non_zero_dual_edg]}
        dic_dual_ver={}  # {key: [dual_ver_coor]}
        map_edg_orig_dic={}  # {(edge): [dual_edge]} #non-zero edges

        for ind, q in enumerate(q_lis):
            den=q[0]  # q is a list of lists
            if abs(den)>self.den_tol:
                edg_dic_new[ind]=edg_dic[ind]
                dic_dual_edg[ind]=dual_edg_lis[ind]
                map_edg_orig_dic[edg_dic[ind]]=dual_edg_lis[ind]
                if edg_dic[ind][0] not in ver_dic_new:
                    ver_dic_new[edg_dic[ind][0]]=ver_dic[edg_dic[ind][0]]
                if edg_dic[ind][1] not in ver_dic_new:
                    ver_dic_new[edg_dic[ind][1]]=ver_dic[edg_dic[ind][1]]
                if dual_edg_lis[ind][0] not in dic_dual_ver:
                    dic_dual_ver[dual_edg_lis[ind][0]]=xyz_s_lis[dual_edg_lis[ind][0]]
                if dual_edg_lis[ind][1] not in dic_dual_ver:
                    dic_dual_ver[dual_edg_lis[ind][1]]=xyz_s_lis[dual_edg_lis[ind][1]]

        self.dic_attr['ver_dic_new']=ver_dic_new  # used in networks_form_duals
        self.dic_attr['edg_dic_new']=edg_dic_new  # used in networks_form_duals
        self.dic_attr['dic_dual_edg']=dic_dual_edg  # used in networks_form_duals
        self.dic_attr['dic_dual_ver']=dic_dual_ver  # used in networks_form_duals
        self.dic_attr['map_edg_orig_dic']=map_edg_orig_dic

    def networks_from_duals(self):
        """
        makes "compas" networks of form and force diagrams (without zero force/cross section elements)
        """
        ver_dic_new=self.dic_attr['ver_dic_new']
        dic_dual_ver=self.dic_attr['dic_dual_ver']
        edg_dic_new=self.dic_attr['edg_dic_new']
        dic_dual_edg=self.dic_attr['dic_dual_edg']
        form_orig_net=hf.make_network(ver_dic_new, edg_dic_new)
        force_orig_net=hf.make_network(dic_dual_ver, dic_dual_edg)

        self.dic_attr['form_orig_net']=form_orig_net
        self.dic_attr['force_orig_net']=force_orig_net

    def process_dual_diagrams(self):
        """
        removes extra vertices and edges from dual diagrams (networks)
        rotates the dual (force) diagram 90 degrees
        calculates the reduired mappings (dictionries) and the updates networks
        saves requires files to draw the processed diagrams in Rhino
        """
        ags_net=self.dic_attr['ags_net']
        form_orig_net=self.dic_attr['form_orig_net']
        force_orig_net=self.dic_attr['force_orig_net']
        map_edg_orig_dic=self.dic_attr['map_edg_orig_dic']
        q_c=self.dic_attr['q_c']  # force_densities, based on dic_attr['edg_dic'] indeces(indeces of original ags_net)
        edg_dic=self.dic_attr['edg_dic']  # the dictionary with original indeces

        # map the original edges to their forces
        old_edg_f_dic={}  # {old_edg:f}
        for ind, edg in edg_dic.items():
            old_q=round(q_c[ind][0], 1)
            old_len=hf.edge_length(ags_net, edg)
            old_edg_f_dic[edg]=(old_q*old_len).item()   # .item() to make it reabale in ironpytho (numpyfloat64>>float)
        
        # update the dual edge mapping (removing repetative vertices of force)
        map_edg_temp_dic=hf.update_dual_mapping_1(force_orig_net, map_edg_orig_dic)

        # update the dual edge mapping
        map_edg_dic, new_edg_f_dic=hf.update_dual_mapping_2(form_orig_net, map_edg_temp_dic, old_edg_f_dic)

        # make a new form_net (without aligned edges)
        form_net=hf.make_new_network(form_orig_net, list(map_edg_dic.keys()))

        # make a new dual (force) network without repetative egdes and vertices
        force_net=hf.make_new_network(force_orig_net, list(map_edg_dic.values()))

        # rotate force_net 90 degrees
        ANG=np.pi/2.0
        force_90_net=hf.rotate_dual(force_net , ANG)

        # dictionary of dual vertices
        dual_ver_dic={}
        for key in force_net.nodes():
            dual_ver_dic[key]=force_net.node_coordinates(key)

        # ### save the data to draw form and force diagrams in Rhino ###
        with open(os.path.join(BASEDIR, 'map_edg_dic.p'), 'wb') as fp:
            pickle.dump(map_edg_dic, fp, protocol=2)
        with open(os.path.join(BASEDIR, 'new_edg_f_dic.p'), 'wb') as fp:
            pickle.dump(new_edg_f_dic, fp, protocol=2)
        with open(os.path.join(BASEDIR, 'dual_ver_dic.p'), 'wb') as fp:
            pickle.dump(dual_ver_dic, fp, protocol=2)          

        self.dic_attr['map_edg_dic']=map_edg_dic
        self.dic_attr['form_net']=form_net
        self.dic_attr['force_net']=force_net
        self.dic_attr['force_90_net']=force_90_net
        self.dic_attr['new_edg_f_dic']=new_edg_f_dic  # {new_edg:f} 

    def constant_stress_fields(self):
        """
        plots and calculates constant stress fields
        """
        new_edg_f_dic=self.dic_attr['new_edg_f_dic']
        map_edg_dic=self.dic_attr['map_edg_dic']
        form_net=self.dic_attr['form_net']
        force_90_net=self.dic_attr['force_90_net']
        form_leaves=form_net.leaves()

        TH=1.0  # [m]  thickness of the coneret element
        SC=1.0/(self.sig_c*TH)  # [m], [(kN/100)/m^2]

        # find convex hulls, anchor lines and memeber CLs for stress field
        hull_lis=[]
        ten_lines_dic={}  # {index:[[(pt1_line1),(pt2_line1)],[(pt1_line2),(pt2_line2)]]}
        ten_ind=0 
        for edg, dual_edg in map_edg_dic.items():
            if edg[0] in form_leaves: # to avoid producing extra anchor points at leaves
                CONST=0
            elif edg[1] in form_leaves:
                CONST=1
            else: 
                CONST=2   
            coor1=form_net.node_coordinates(edg[0])
            coor2=form_net.node_coordinates(edg[1])
            dual_coor1=force_90_net.node_coordinates(dual_edg[0])
            dual_coor2=force_90_net.node_coordinates(dual_edg[1])
            sc_pts_lis=scale_points([dual_coor1, dual_coor2], SC)
            if new_edg_f_dic[edg]>self.den_tol:  # tension
                line_pts_lis=hf.sf_cl_anchor_lines(sc_pts_lis, [coor1, coor2], CONST)
                ten_lines_dic[ten_ind]=line_pts_lis
                ten_ind+=1
            elif new_edg_f_dic[edg]<-self.den_tol:  # compression
                hull=hf.minkowski_sum(sc_pts_lis, [coor1, coor2])
                hull_lis.append(hull)

        hf.plot_stress_fields(hull_lis, ten_lines_dic)
        hf.sf_rhino_inputs(hull_lis, ten_lines_dic)

"################################################ INPUTS #################################################"
# mesh_cant (cantilever beam-Fig.6)
dic_sup={21:[1, 1], 7:[1, 1]}
dic_load={20:[0.0, -100.0]}

# mesh_abut (abutment-Fig.7)
# dic_sup={48:[1, 1], 4:[0, 1]}
# dic_load={55:[-50.0, 0.0]}

# mesh_ssb (simply-supported beam-Fig.9)
# dic_sup={0:[1, 1], 12:[0, 1]}  
# dic_load={56:[0.0, -200.0], 60:[0.0, -200.0]} 

# mesh_step (stepped beam-Fig.9)
# dic_sup={52:[1, 1], 32:[1, 0]}
# dic_load={31:[400.0, 0.0], 7:[-400.0, 0.0]}

# mesh_joint (frame joint-Fig.9)
# dic_sup={30:[1, 1], 0:[0, 1]}
# dic_load={40:[400, 0.0], 43:[-400.0, 00.0]}

# mesh_dap (depped-end beam-Fig.9)
# dic_sup={12:[0, 1], 0:[1, 1]}
# dic_load={9:[-250.0, 0.0], 6:[250.0, 0.0], 7:[-75.0, -75.0]} 

# mesh_hhp (hammerhead pier-Fig.9)
# dic_sup={33:[1, 1], 37:[1, 1]}
# dic_load={23:[0.0, -80.0], 27:[0.0, -30.0], 31:[0.0, -10.0]}

# mesh_wall (wall with opening-Fig.9)
# dic_sup={3:[1, 1], 67:[0, 1]} 
# dic_load={61:[0.0, -400.0]}

"################################################ RUN ###################################################"

#************ get the initial meshed concrete block and create ground truss ************
AO=AGS_OPT_2D(dic_sup, dic_load)
# !!! ENTER THE CORRECT DIRECTORY BELOW: !!!
AO.create_mesh(r'C:\Users\...\MeshObjects\mesh_cant.obj')
# hf.plot_mesh(AO.dic_attr['mesh'])  #!uncomment to see the input mesh
AO.struct_domain(AO.dic_attr['mesh'])
AO.generate_ground_truss()
AO.network_from_ground_truss()
# hf.plot_network(AO.dic_attr['gt_net']) #!uncomment to see the ground truss

# ************ add vertices at intersections and add load/sup edges to gt_net  ************
ags_net=AO.post_processing(False, False, True, True, True, AO.dic_attr['gt_net'])
AO.dic_attr['ags_net']=ags_net
# hf.plot_network(ags_net)

# ************ create a new GT with removed load/sup edges, get new GT properties with mid-nodes ************
gtopt_net=AO.post_processing(False, True, False, False, False, ags_net)
AO.dic_attr['gtopt_net']=gtopt_net
# hf.plot_network(gtopt_net)

# *********** inputs for AGS ************
with open(os.path.join(BASEDIR, 'ver_dic_GT.p'), 'rb') as fp:
    AO.dic_attr['ver_dic']=pickle.load(fp)
with open(os.path.join(BASEDIR, 'edg_dic_GT.p'), 'rb') as fp:
    AO.dic_attr['edg_dic']=pickle.load(fp)

# ************ initialize AGS >> OPT ************
AO.inc_mat()
AO.inc_s_mat()
AO.coor_dif()
AO.equil_mat()
AO.fixed_dof()
AO.free_dof()
AO.load_vec()
AO.lin_prog_inputs()
AO.lin_prog()
AO.leaf_edg_dof_mapping()
AO.make_q_complete()
AO.dual_coor()
AO.plot_dual_diagrams()

# ************ post-process the dual diagrams ************
AO.update_duals_info()
AO.networks_from_duals()
AO.process_dual_diagrams()
AO.constant_stress_fields()

#************ show plots of form, force, and "unadjusted" stress fields ************
plt.show()

#************ print force magnitudes ************
# !!! CLOSE PREVIOUS PLOTS !!!
# hf.plot_network(AO.dic_attr['form_net']) # to see the vertex indices
# print (AO.dic_attr['new_edg_f_dic']) # negative=compression, positive=tension
