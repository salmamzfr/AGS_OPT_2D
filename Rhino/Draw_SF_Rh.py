'''
Created on 20.06.2018, updated on 23.01.2020

@author: Salma Mozaffari, ETH Zurich

Before running the code:
    update the directories below "import files" according to your directory.
'''
import rhinoscriptsyntax as rs
import cPickle as pickle

if not rs.IsLayer("anchor"):
    rs.AddLayer("anchor", (0, 0, 0))
if not rs.IsLayer("tie"):
    rs.AddLayer("tie", (255, 0, 0))
if not rs.IsLayer("sf"):
    rs.AddLayer("sf", (0, 0, 225))

SC=(1.0, 1.0, 1.0)
ORIG=(0.0, 0.0, 0.0)

# import files !!! ENTER THE CORRECT DIRECTORY BELOW: !!!
with open(r'C:\Users\...\AGS_OPT_2D\Source\hull_dic.p', 'rb') as fp:
    hull_dic=pickle.load(fp)
with open(r'C:\Users\...\AGS_OPT_2D\Source\ten_lines_dic.p', 'rb') as fp:
    ten_lines_dic=pickle.load(fp)

# draw compresson stress fields
for seg_lis in hull_dic.values():
    for seg in seg_lis:
        pt_1=(seg[0][0], seg[0][1], 0.0)
        pt_2=(seg[1][0], seg[1][1], 0.0)
        ln=rs.AddLine(pt_1,pt_2)
        object_id = rs.FirstObject()
        rs.ScaleObject(object_id, ORIG, SC, copy=False)
        rs.ObjectLayer(rs.FirstObject(),"sf")

# draw tension ties and anchoring regions
for lines in ten_lines_dic.values():
    line_1=lines[0]
    pt_1=(line_1[0][0], line_1[0][1], 0.0)
    pt_2=(line_1[1][0], line_1[1][1], 0.0)
    ln=rs.AddLine(pt_1,pt_2)
    object_id = rs.FirstObject()
    rs.ScaleObject(object_id, ORIG, SC, copy=False)
    rs.ObjectLayer(rs.FirstObject(),"tie")
    lines_23 =lines[1:]
    for line in lines_23:
        pt_1=(line[0][0], line[0][1], 0.0)
        pt_2=(line[1][0], line[1][1], 0.0)
        ln=rs.AddLine(pt_1,pt_2)
        object_id = rs.FirstObject()
        rs.ScaleObject(object_id, ORIG, SC, copy=False)
        rs.ObjectLayer(rs.FirstObject(),"anchor")