import numpy as np
import os
import math
import random
from numpy import genfromtxt
import numpy as npimport 
from numpy import genfromtxt
#import networkx as nx
#import matplotlib.pyplot as plt  # export DISPLAY=localst:0.0
#from pysmiles import read_smiles
#import spektral
#from spektral import chem
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PyMol
# could use novel_score from spektral
#from pysmiles import read_smiles
from numpy import savetxt

# 1's are used to indicate rings. must have following first atom and following last atom that completes the ring
"""
G = nx.Graph()
mol_with_H = read_smiles(m, explicit_hydrogen=True)
nx.draw(mol_with_H)
mol_with_H = read_smiles('OC(=O)C(CNC)C(C=O)', explicit_hydrogen=True)
nx.draw(mol_with_H)
plt.show()
"""
# 'OC(=O)O'

element_dict = {
    "H"  : 1,
    "He" : 2,
    "Li" : 3,
    "Be" : 4,
    "B"  : 5,
    "C"  : 6,
    "N"  : 7,
    "O"  : 8,
    "F"  : 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P" : 15,
    "S" : 16,
    "Cl": 17,
    "Ar": 18,
    "K" : 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V" : 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y" : 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I" : 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W" : 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U" : 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm":100,
    "Md":101,
    "No":102,
    "Lr":103,
    "Rf":104,
    "Db":105,
    "Sg":106,
    "Bh":107,
    "Hs":108,
    "Mt":109,
    "Ds":110,
    "Rg":111,
    "Cn":112,
    "Nh":113,
    "Fl":114,
    "Mc":115,
    "Lv":116,
    "Ts":117,
    "Og":118

}

def distance(geom1, geom2):
    """ Three dimensional distance formula for evaluating the distance betweeen molecules """
    return math.sqrt( (geom1[1] - geom2[1]) **2 + (geom1[2] - geom2[2])**2 + ( geom1[3] - geom2[3])**2)

def convertTuple(tup):
    str =  ''.join(tup)
    return str

def bond_lengths_no_f(geom):

    lines = []
    for i in range(len(geom[:,0])):
        for j in range(len(geom[:,0])):
                
            distances = round(distance(geom[i,:], geom[j,:]), 3)

            if j > i and distances < 3.0 :
        
                line = [i + 1, j + 1, distances]

                lines.append(line)

            if j > i and distances < 3.0 and i == len(geom[:,0]):
        
                line = [i + 1, j + 1, distances]

                lines.append(line)
    array = np.array(lines)

    #with open ('sec_1.txt', 'w') as fp: 
    #    fp.write(lines)
    return array
           
    #with open ('bonds.txt', 'w') as fp: 
    #    fp.write(lines)

    return lines

def bond_angles_no_f(geom):
    length = len(geom[:,0])
    angles = []
    for i in range(length - 2):

            ab = distance(geom[i,:], geom[i+1,:])
            ac = distance(geom[i,:], geom[i+2,:])
            bc = distance(geom[i+1,:], geom[i+2,:])

            angle = math.degrees(math.acos( ( ab **2 + ac **2 - bc **2 )/(2*ab*ac)))

            angle = round(angle, 3)
            if i < length:
                ang = str(i + 1)+ " ", str(i +2)+ " ", str(i +3)+ " ", "=" +str(angle), " B\n", str(i + 1)+ " ", str(i +2)+ " ", str(i +3)+ " F\n"
                ang = convertTuple(ang)
                angles.append(ang)
            if i == length:
                ang = str(i + 1)+ " ", str(i +2)+ " ", str(i +3)+ " ", "=" +str(angle), " B\n", str(i + 1)+ " ", str(i +2)+ " ", str(i +3)+ " F"
                ang = convertTuple(ang)
                angles.append(ang)
    angles = ''.join(angles)            
    #with open ('angles.txt', 'w') as fp: 
    #
    #     fp.write(angles)
    return angles

def dihedral_no_f(geom):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    length = len(geom[:,0])
    di_angles = []

    for i in range(length - 3):
        p0 = geom[i,1:] # in the form np.array[x,y,z]
        p1 = geom[i+1,1:]
        p2 = geom[i+2,1:]
        p3 = geom[i+3,1:]

        b0 = -1.0*(p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2

        # normalize b1 so that it does not influence magnitude of vector
        # rejections that come next
        b1 /= np.linalg.norm(b1)

        # vector rejections
        # v = projection of b0 onto plane perpendicular to b1
        #   = b0 minus component that aligns with b1
        # w = projection of b2 onto plane perpendicular to b1
        #   = b2 minus component that aligns with b1
        v = b0 - np.dot(b0, b1)*b1
        w = b2 - np.dot(b2, b1)*b1

        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        di_ang = np.degrees(np.arctan2(y, x))
        di_angle = round(di_ang, 3)
        if i < length:
            di_ang = (str(i + 1)+ " ", str(i +2)+ " ", str(i +3)+ " ", str(i +4)+ " ", "=" +str(di_angle), " B\n", 
                str(i + 1)+ " ", str(i +2)+ " ", str(i +3) + " " + str(i +4)+ " "," F\n")
            di_ang = convertTuple(di_ang)
            di_angles.append(di_ang)
        if i == length:
            di_ang = (str(i + 1)+ " ", str(i +2)+ " ", str(i +3)+ " ", str(i +4)+ " ", "=" +str(di_angle), " B\n", 
                str(i + 1)+ " ", str(i +2)+ " ", str(i +3) + " " + str(i +4)+ " "," F")
            di_ang = convertTuple(di_ang)
            di_angles.append(di_ang)
    
    di_angles = ''.join(di_angles)            
    #with open ('dihedral.txt', 'w') as fp: 
    #    fp.write(di_angles)
    return di_angles

def clean_carts(num_u_confs, group=0):

    """ This will replace the numerical forms of the elements as their letters numbered in order """
    for m in range(num_u_confs):
        filename = "carts_" + str(group*num_u_confs + m +1) + ".txt"
        f = open(filename,'r')
        a = ["118.0 ",
            "117.0 ", 
            "116.0 ", 
            "115.0 ", 
            "114.0 ", 
            "113.0 ", 
            "112.0 ", 
            "111.0 ", 
            "110.0 ", 
            "109.0 ", 
            "108.0 ", 
            "107.0 ", 
            "106.0 ", 
            "105.0 ", 
            "104.0 ", 
            "103.0 ", 
            "102.0 ", 
            "101.0 ", 
            "100.0 ", 
            "99.0 ", 
            "98.0 ", 
            "97.0 ", 
            "96.0 ", 
            "95.0 ", 
            "94.0 ", 
            "93.0 ", 
            "92.0 ", 
            "91.0 ", 
            "90.0 ", 
            "89.0 ", 
            "88.0 ", 
            "87.0 ", 
            "86.0 ", 
            "85.0 ", 
            "84.0 ", 
            "83.0 ", 
            "82.0 ", 
            "81.0 ", 
            "80.0 ", 
            "79.0 ", 
            "78.0 ", 
            "77.0 ", 
            "76.0 ", 
            "75.0 ", 
            "74.0 ", 
            "73.0 ", 
            "72.0 ", 
            "71.0 ", 
            "70.0 ", 
            "69.0 ", 
            "68.0 ", 
            "67.0 ", 
            "66.0 ", 
            "65.0 ", 
            "64.0 ", 
            "63.0 ", 
            "62.0 ", 
            "61.0 ", 
            "60.0 ", 
            "59.0 ", 
            "58.0 ", 
            "57.0 ", 
            "56.0 ", 
            "55.0 ", 
            "54.0 ", 
            "53.0 ", 
            "52.0 ", 
            "51.0 ", 
            "50.0 ", 
            "49.0 ", 
            "48.0 ", 
            "47.0 ", 
            "46.0 ", 
            "45.0 ", 
            "44.0 ", 
            "43.0 ", 
            "42.0 ", 
            "41.0 ", 
            "40.0 ", 
            "39.0 ", 
            "38.0 ", 
            "37.0 ", 
            "36.0 ", 
            "35.0 ", 
            "34.0 ", 
            "33.0 ", 
            "32.0 ", 
            "31.0 ", 
            "30.0 ", 
            "29.0 ", 
            "28.0 ", 
            "27.0 ", 
            "26.0 ", 
            "25.0 ", 
            "24.0 ", 
            "23.0 ", 
            "22.0 ", 
            "21.0 ", 
            "20.0 ", 
            "19.0 ", 
            "18.0 ", 
            "17.0 ", 
            "16.0 ", 
            "15.0 ", 
            "14.0 ", 
            "13.0 ", 
            "12.0 ", 
            "11.0 ", 
            "10.0 ", 
            "9.0 ",
            "8.0 ",
            "7.0 ",
            "6.0 ",
            "5.0 ",
            "4.0 ",
            "3.0 ",
            "2.0 ",
            "1.0 ", 


            ]
        table = {"1.0 "  :"H"  , 
                "2.0 "  :"He" ,
                "3.0 "  :"Li" ,
                "4.0 "  :"Be" ,
                "5.0 "  :"B"  ,
                "6.0 "  :"C"  ,
                "7.0 "  :"N"  ,
                "8.0 "  :"O"  ,
                "9.0 "  :"F"  ,
                "10.0 " :"Ne" ,
                "11.0 " :"Na" ,
                "12.0 " :"Mg" ,
                "13.0 " :"Al" ,
                "14.0 " :"Si" ,
                "15.0 " :"P"  ,
                "16.0 " :"S"  ,
                "17.0 " :"Cl" ,
                "18.0 " :"Ar" ,
                "19.0 " :"K"  ,
                "20.0 " :"Ca" ,
                "21.0 " :"Sc" ,
                "22.0 " :"Ti" ,
                "23.0 " :"V"  ,
                "24.0 " :"Cr" ,
                "25.0 " :"Mn" ,
                "26.0 " :"Fe" ,
                "27.0 " :"Co" ,
                "28.0 " :"Ni" ,
                "29.0 " :"Cu" ,
                "30.0 " :"Zn" ,
                "31.0 " :"Ga" ,
                "32.0 " :"Ge" ,
                "33.0 " :"As" ,
                "34.0 " :"Se" ,
                "35.0 " :"Br" ,
                "36.0 " :"Kr" ,
                "37.0 " :"Rb" ,
                "38.0 " :"Sr" ,
                "39.0 " :"Y"  ,
                "40.0 " :"Zr" ,
                "41.0 " :"Nb" ,
                "42.0 " :"Mo" ,
                "43.0 " :"Tc" ,
                "44.0 " :"Ru" ,
                "45.0 " :"Rh" ,
                "46.0 " :"Pd" ,
                "47.0 " :"Ag" ,
                "48.0 " :"Cd" ,
                "49.0 " :"In" ,
                "50.0 " :"Sn" ,
                "51.0 " :"Sb" , 
                "52.0 " :"Te" , 
                "53.0 " :"I"  , 
                "54.0 " :"Xe" , 
                "55.0 " :"Cs" , 
                "56.0 " :"Ba" , 
                "57.0 " :"La" , 
                "58.0 " :"Ce" , 
                "59.0 " :"Pr" , 
                "60.0 " :"Nd" , 
                "61.0 " :"Pm" , 
                "62.0 " :"Sm" , 
                "63.0 " :"Eu" , 
                "64.0 " :"Gd" , 
                "65.0 " :"Tb" , 
                "66.0 " :"Dy" , 
                "67.0 " :"Ho" , 
                "68.0 " :"Er" , 
                "69.0 " :"Tm" , 
                "70.0 " :"Yb" ,  
                "71.0 " :"Lu" ,  
                "72.0 " :"Hf" ,  
                "73.0 " :"Ta" ,  
                "74.0 " :"W"  ,  
                "75.0 " :"Re" ,  
                "76.0 " :"Os" ,  
                "77.0 " :"Ir" ,  
                "78.0 " :"Pt" ,  
                "79.0 " :"Au" ,  
                "80.0 " :"Hg" ,  
                "81.0 " :"Tl" ,  
                "82.0 " :"Pb" ,  
                "83.0 " :"Bi" ,  
                "84.0 " :"Po" ,  
                "85.0 " :"At" ,  
                "86.0 " :"Rn" ,  
                "87.0 " :"Fr" ,   
                "88.0 " :"Ra" ,   
                "89.0 " :"Ac" ,   
                "90.0 " :"Th" ,   
                "91.0 " :"Pa" ,   
                "92.0 " :"U"  ,   
                "93.0 " :"Np" , 
                "94.0 " :"Pu" , 
                "95.0 " :"Am" , 
                "96.0 " :"Cm" , 
                "97.0 " :"Bk" , 
                "98.0 " :"Cf" , 
                "99.0 " :"Es" , 
                "100.0 " :"Fm", 
                "101.0 " :"Md", 
                "102.0 " :"No", 
                "103.0 " :"Lr", 
                "104.0 " :"Rf", 
                "105.0 " :"Db", 
                "106.0 " :"Sg", 
                "107.0 " :"Bh", 
                "108.0 " :"Hs", 
                "109.0 " :"Mt",  
                "110.0 " :"Ds",  
                "111.0 " :"Rg",  
                "112.0 " :"Cn",  
                "113.0 " :"Nh",  
                "114.0 " :"Fl",  
                "115.0 " :"Mc",  
                "116.0 " :"Lv",  
                "117.0 " :"Ts",  
                "118.0 " :"Og"   

        } # this method is finding the single digit values before the double digit for tranlation

        lst = []
        cnt2 = 0
        for line in f:
            cnt2 += 1
            for word in a: 
                if word in line:
                    convert_wrd = table[word]
                    line = line.replace(word, convert_wrd + str(cnt2) + " ")
                    
            lst.append(line)
        f.close()
        f = open(filename,'w')
        for line in lst:
            f.write(line)
        f.close()


def make_conf_input_dir(num_u_confs, group=0):
    """ Combines the geometry output and the constrained output. Then makes the .com and .pbs files in a subdirectory """
    for m in range(num_u_confs):
        data = ""
        filename = "carts_" + str(group*num_u_confs + m +1) + ".txt"
        with open(filename) as fp: 
            data = fp.read() 
        

        data += "\n\n"
        charges = "0 1"

        new_dir = "geom" + str(group + m+1)
        os.mkdir(new_dir)
        with open (new_dir + '/mex.com', 'w') as fp: 
            fp.write("%mem=1600mb\n")
            fp.write("%nprocs=4\n")
            fp.write("#N wB97XD/6-31G(d) opt=ModRedundant FREQ\n")
            fp.write("\n")
            fp.write("Name ModRedundant\n")
            fp.write("\n")
            fp.write(charges + "\n")
            fp.write(data) 

        with open (new_dir + '/mex.pbs', 'w') as fp: 
            fp.write("#!/bin/sh\n")
            fp.write("#PBS -N mex\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m abe\n#PBS -l")
            fp.write("mem=15gb\n")
            fp.write("#PBS -l nodes=1:ppn=4\n#PBS -q gpu\n\nscrdir=/tmp/$USER.$PBS_JOBID\n\n")
            fp.write("mkdir -p $scrdir\nexport GAUSS_SCRDIR=$scrdir\nexport OMP_NUM_THREADS=1\n\n")
            fp.write("""echo "exec_host = $HOSTNAME"\n\nif [[ $HOSTNAME =~ cn([0-9]{3}) ]];\n""")
            fp.write("then\n")
            fp.write("  nodenum=${BASH_REMATCH[1]};\n  nodenum=$((10#$nodenum));\n  echo $nodenum\n\n")
            fp.write("""  if (( $nodenum <= 29 ))\n  then\n    echo "Using AVX version";\n""")
            fp.write("    export g16root=/usr/local/apps/gaussian/g16-b01-avx/\n  elif (( $nodenum > 29 ))\n")
            fp.write("""  then\n    echo "Using AVX2 version";\n    export g16root=/usr/local/apps/gaussian/g16-b01-avx2/\n  else\n""")
            fp.write("""    echo "Unexpected condition!"\n    exit 1;\n  fi\nelse\n""")
            fp.write("""  echo "Not on a compute node!"\n  exit 1;\nfi\n\n""")
            fp.write("cd $PBS_O_WORKDIR\n. $g16root/g16/bsd/g16.profile\ng16 mex.com mex.out\n\nrm -r $scrdir\n")

        #print("\ngeom" + str(dir_name_number) + "\n")
        #os.chdir("geom" + str(dir_name_number))
        #os.system("qsub mex.pbs")
        #os.chdir("..")

def build_confs(smiles, conf_num=10, element_dict=element_dict):

    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m, addCoords=True)   # uncomment for hydrogen filling


    num_atoms = m.GetNumAtoms()
    AllChem.EmbedMultipleConfs(m,conf_num)
    rmslist = []
    AllChem.AlignMolConformers(m, RMSlist=rmslist)
    #print(rmslist) rms values from first conformer to the others. AllChem.GetConformerRMS(m2, 1, 9, prealigned=True) for others

    for i in range(m.GetNumConformers()):
        AllChem.UFFOptimizeMolecule(m,confId=i)
        
        #v.ShowMol(m,confId=i,name='conf-%d'%i,showOnly=False)
    w = Chem.SDWriter('confs_sdf.sdf')
    for i in range(m.GetNumConformers()):
        w.write(m,confId=i)
        

    w.flush()

    f=open('confs_sdf.sdf','r')
    lines = f.readlines()
    f.close()


    baseline_in = False

    breaker = "$$$$"
    rm_lines= []
    lines = []

    with open('confs_sdf.sdf') as search:
        for num, line in enumerate(search,1):
            lines.append(line)
            if breaker in line:
                rm_lines.append(num)
                baseline_in = True

    with open('confs_sdf.sdf') as f:
        lines = f.read().splitlines()


    new_ls = []
    for i in lines:
        cart_bonds = ' '.join(i.split())
        new_ls.append(cart_bonds)


    sdf_dict = {}
    for k in range(1, conf_num+1):
        if k == 1:
            sdf_dict['sec_{0}'.format(k)] = (new_ls[4:rm_lines[0]-2])
        elif k > 1 and k < conf_num +1:
            sdf_dict['sec_{0}'.format(k)] = new_ls[rm_lines[k-2] +4:rm_lines[k-1] - 2]
    return sdf_dict, num_atoms, conf_num



def sdf_dict_cart(sdf_dict, num_atoms, element_dict, conf_num, group=0):
    bond_length = {}
    cnt = 1

    for key in sdf_dict:
        ls = []
        cur_sec = sdf_dict[key]
        for i in range(len(cur_sec)):
            ls.append(cur_sec[i].split(' '))
        
        
        ar_c = np.zeros((num_atoms, 4))
        ar_i = np.zeros((len(ls[num_atoms:]), 4))
        
        for num, j in enumerate(ls[:num_atoms]):
            
            ls[num][3] = element_dict[j[3]]    
        
        for i in range(num_atoms):


            ar_c[i,0] = ls[i][3]
            ar_c[i,1] = ls[i][0]
            ar_c[i,2] = ls[i][1]
            ar_c[i,3] = ls[i][2]

        if key == 'sec_1':
            arching = ar_c
            #print(ar_c) S still here
        else:
            arching = np.concatenate((arching, ar_c), axis=0)
    

        bond_length['sec_{0}'.format(cnt)] = bond_lengths_no_f(ar_c)
        cnt +=1
    
    out_file = "carts.txt"

    np.savetxt(out_file, arching,
    fmt="%s")

     
    check_list = [] 
    final_dict = {}
    final_dict['sec_1'] = sdf_dict['sec_1']

    t = 0
    for num_1, k_1 in enumerate(bond_length):
        for num_2, k_2 in enumerate(bond_length):
            if k_1 != k_2 and num_2 > num_1: # creates a series of [n, n+1 : (-n)]
                in_check = []
                #print(k_1, k_2)
                for n_1, i in enumerate(bond_length[k_1]):
                    for n_2, j in enumerate(bond_length[k_2]):
                        t += 1
                        if n_1 == n_2:
                            if abs( round(i[2], 2) - round(j[2], 2)) < 0.08:
                                in_check.append(True)
                            else: 
                                in_check.append(False)
                if np.all(in_check) == True:
                    check_list.append((num_1, num_2, 1))
                else:
                    check_list.append((num_1, num_2, 0)) 
                    
    tripped = []
    delete = []

    for i in check_list:
        pair = '{}{}'.format(i[0], i[1])
        if i[2] == 1:
            tripped.append(pair)
            delete.append(i[1]) # shows which are duplicates
    delete = list(set(delete))
    #print(delete)

    del_cnt = 0
    for k in range(conf_num):
        #print(k)
        if k in delete:
            #print(k, "deleted")
            del_cnt += 1
        else:
            if k == 0:
                final = arching[:num_atoms,:]
                #print(final)
            elif k < conf_num -1:
                final = np.concatenate((final, arching[k*num_atoms : k*num_atoms + num_atoms, :]))

            elif k == conf_num:
                final = np.concatenate((final, arching[k*num_atoms :, :]))


    
    num_u_confs = conf_num - del_cnt
    
    out_file = "carts_no_confs.txt"

    np.savetxt(out_file, final,
    fmt="%s")

    for m in range(num_u_confs):
        out_file = "carts_" + str(group*num_u_confs + m +1) + ".txt"
        
        piece = final[m*num_atoms: (m+1)*num_atoms]
        #if m == 1:
            #print(piece) 
        np.savetxt(out_file, piece,
        fmt="%s")

    clean_carts(num_u_confs, group)
    #make_conf_input_dir(num_u_confs)

    #print(num_u_confs) 
    return num_u_confs


# single - (or blank), double = , triple #

#smiles = 'OC(=O)C(C=O)(C1NC)C(NC4CCCCCCC3CCCC(CCCCCC(CCC(CCC4CC)C)CCCCCCCC)CCC2)(CCC(C=O)NO2)COOC(CCCC3C)OOO(C1=O)'
#smiles = 'C1CCCCC1N(C)C(=O)C(=CI)C(O)C=C'
#smiles = 'C1(C(O)(C)(C))=CC=CC=C1CCC(SCC2CC2CC(=O)O)C3=CC=CC(=C3)C=CC4=CC=C(C5=N4)C=CC(Cl)=C5' #Montelukast
smiles = 'C1=CC=CC=C1C(=O)NC(C2=CC=CC=C2)C(O)C(=O)OC(C(C)-C3(C(OC(=O)C))C(=O)C4(C)C(O)CC(OC5)C5(OC(=O)C)C4C7(OC(=O)C6=CC=CC=C6))CC7(O)(C3(C)(C))' # Paclitaxel
#smiles = 'C1CCCC1CSCCCCC'
#smiles = 'OC(=O)O'
#smiles = 'NC(B)(Cl)O'

""" def randomSmiles(m1):
    m1.SetProp("_canonicalRankingNumbers", "True")
    idxs = list(range(0,m1.GetNumAtoms()))
    random.shuffle(idxs)
    for i,v in enumerate(idxs):
        m1.GetAtomWithIdx(i).SetProp("_canonicalRankingNumber", str(v))
    return Chem.MolToSmiles(m1)

def molecular_formula_input(smiles):

    m1 = Chem.MolFromSmiles(smiles)
    print(m1)
    s = set()
    for i in range(1000):
        smiles = randomSmiles(m1)
        s.add(smiles)
    s = list(s)
    return s



const_iso = molecular_formula_input(smiles) """

input_mo_form = "N1C1B1O1"

def molecular_formula_input(molecular_formula):
    input_mo_form = "H2C1O2"
    input_mo_form = "N1C1B1O1"

    def split(word):
        return [char for char in word]


    def atom_list(split_in):
        cnt = 1
        new = []
        net = []
        for i in split_in:
            if cnt%2 == 2%2:
                new.append(int(i ))
            else:
                net.append(i)
            cnt += 1

        nex = ""
        for i, j in zip(net, new):
            nex = nex + i*j
        
        nex = split(nex)
        neb = []
        for i, j in zip(nex, range(len(nex))):
            neb.append(i)

        return neb

    split_in = split(input_mo_form)

    ele_s = atom_list(split_in)
    print( atom_list(split_in))

    #for el in atom_list(split_in):
    #   print(el, type(el))

    smile_time = []
    smiles = set()
    middle = 0

    for i in range(1):
        neb = ele_s[:]
        smile_time = []
        for j in range(len(neb)):
            print(j)
            rand_i = np.random.choice(neb)
            paren = np.random.choice([0,1,2,3])
            if paren == 1 and j != 0:
                smile_time.append("(" + rand_i)
                neb.remove(rand_i)
                middle += 1

            elif j==len(neb) - 1 and middle > 0:
                print('yes')
                for k in range(middle +1):
                    smile_time.append(rand_i + ")")
                    neb.remove(rand_i)
            
            elif paren == 1 and middle > 0:
                smile_time.append( rand_i + ")")
                neb.remove(rand_i)
                middle += -1
            

            else:
                smile_time.append(rand_i)
                neb.remove(rand_i)
            print(smile_time)
            # could have running total of plus one for left parentheses and minus one for right. If at the end the value is not zero, then need to add parentheses
            

        smile_time = ''.join(smile_time)
        smiles.add(smile_time) # produces n! structures, need to find ways to cut down generations
    smiles = list(smiles)
    #smiles is a set of all permutations
    print(smiles)
    return smiles

#const_iso = molecular_formula_input(input_mo_form) # molecular formula

const_iso = [smiles] # smiles string connectivity defined

for num, i in enumerate(const_iso):
    group = num
    print(i)
    sdf_dict, num_atoms, conf_num = build_confs(smiles=i, conf_num=1, element_dict=element_dict)
    sdf_dict_cart(sdf_dict, num_atoms, element_dict, conf_num, group=num)
    # need to build constrained molecular formula generator




#sdf_dict, num_atoms, conf_num = build_confs(smiles=smiles, conf_num=10, element_dict=element_dict)
#sdf_dict_cart(sdf_dict, num_atoms, element_dict, conf_num)


