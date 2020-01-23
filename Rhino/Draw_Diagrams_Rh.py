'''
Created on 09.05.2018, updated on 22.01.2020

@author: Salma Mozaffari, ETH Zurich

Before running the code:
    update the directories below "import files" according to your directory.
'''
import rhinoscriptsyntax as rs
import cPickle as pickle

if not rs.IsLayer("comp"):
    rs.AddLayer("comp", (0, 0, 255))
if not rs.IsLayer("ten"):
    rs.AddLayer("ten", (255, 0, 0))

# import files !!! ENTER THE CORRECT DIRECTORY BELOW: !!!
with open(r'C:\Users\...\AGS_OPT_2D\Source\ver_dic_GT.p', 'rb') as fp:
    ver_dic=pickle.load(fp)
with open(r'C:\Users\...\AGS_OPT_2D\Source\new_edg_f_dic.p', 'rb') as fp:
    edg_f_dic=pickle.load(fp)
with open(r'C:\Users\...\AGS_OPT_2D\Source\map_edg_dic.p', 'rb') as fp:
    map_edg_dic=pickle.load(fp)
with open(r'C:\Users\msalma\...\AGS_OPT_2D\Source\dual_ver_dic.p', 'rb') as fp:
    dual_ver_dic=pickle.load(fp)


# draw form and force diagrams (scales adjustable!)
for edg, f in edg_f_dic.items():
    dl_edg=map_edg_dic[edg]
    pt_1, pt_2=[ver_dic[edg[0]][0], ver_dic[edg[0]][1], 0.0], [ver_dic[edg[1]][0], ver_dic[edg[1]][1], 0.0]
    dl_pt_1, dl_pt_2=[dual_ver_dic[dl_edg[0]][0], dual_ver_dic[dl_edg[0]][1], 0.0], [dual_ver_dic[dl_edg[1]][0], dual_ver_dic[dl_edg[1]][1], 0.0]
    if f>=0.01:  # tension
        line_1=rs.AddLine(pt_1, pt_2)
        line_2=rs.AddLine(dl_pt_1, dl_pt_2)
        rs.ObjectLayer([line_1, line_2], "ten")
    elif f<=-0.01:  # compression
        line_1=rs.AddLine(pt_1, pt_2)
        line_2=rs.AddLine(dl_pt_1, dl_pt_2)
        rs.ObjectLayer([line_1, line_2], "comp")

# draw "unadjusted" stress fields

