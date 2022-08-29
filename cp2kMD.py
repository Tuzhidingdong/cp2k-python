class cp2kMD():
    def __init__(self,dire):

        self.dire = dire

        figsizeinch=(3.3464566929,2.5098425197)
        self.figsizeinch=figsizeinch
        self.cm_inch=2.539999918
        self.inch_cm=0.394
        self.dire = dire
        figsizeinch=(3.3464566929,2.5098425197)
        self.figsizeinch=figsizeinch
        threefig=(3.3464566929*0.65,2.5098425197*0.65)
        self.threefig=threefig
        insert=(3.3464566929*0.45,2.5098425197*0.45)
        self.insert=insert
        figsizeinchdouble=(3.3464566929*1.5,2.5098425197)
        self.figsizeinchdouble=figsizeinchdouble
        figsizeneber=(3.3464566929*0.7,2.5098425197)
        self.figsizeneber=figsizeneber
        figsizeneber=(3.3464566929*0.7,2.5098425197)
        self.figsizeneber=figsizeneber        
    def POSTAIMD(self,filename,rdf_fir,rdf_sec,):

        from ase.io import read, write
        import os
        import MDAnalysis
        import MDAnalysis.analysis.rdf
        import MDAnalysis.analysis.rms
        import matplotlib.pyplot as plt
        import pandas as pd


        os.chdir(self.dire)
        # write('traj.pdb', read(filename))



        u = MDAnalysis.Universe(filename, permissive=True)
        u.dimensions = [73, 73, 73, 90, 90, 90]
        u.atoms.write("traj.pdb",frames='all')
        u.atoms.write("crystal.xyz")
        g1= u.select_atoms(rdf_fir)
        g2= u.select_atoms(rdf_sec)
        rdf = MDAnalysis.analysis.rdf.InterRDF(g1,g2,nbins=75, range=(0.0, min(u.dimensions[:3])/2.0))
                  
        rdf.run()
        fig = plt.figure(figsize=(5,4))
        ax = fig.add_subplot(111)
        ax.plot(rdf.bins, rdf.rdf, 'k-',  label="rdf")
        ax.legend(loc="best")
        ax.set_xlabel(r"Distance ($\AA$)")
        ax.set_ylabel(r"RDF")
        fig.savefig("RDF_all.png")





        ref = MDAnalysis.Universe(filename, permissive=True) 
        ref.trajectory[0]
        R = MDAnalysis.analysis.rms.RMSD(u, ref,
           select="all",         # superimpose on whole backbone of all atoms # align based on all atoms
           groupselections=[rdf_fir,rdf_sec],
           filename="rmsd_all.dat",center=True)#,   # CORE
        timestep=0.0005  #0.5fs from fs to ps as Reader has no dt information, set to 1.0 ps          
        R.run()
        rmsd = R.rmsd.T   # transpose makes it easier for plotting
        time = rmsd[1]*timestep

        da = {'time':time,'all':rmsd[2],rdf_fir:rmsd[3],rdf_sec:rmsd[4]}
        df=pd.DataFrame(da)#构造原始数据文件
        df.to_excel("RMSD.xlsx")#生成Excel文件，并存到指定文件路径下




        fig = plt.figure(figsize=(5,4))
        ax = fig.add_subplot(111)
        ax.plot(time, rmsd[2], 'k-',  label="all")
        ax.plot(time, rmsd[3], 'r--', label=rdf_fir)
        ax.plot(time, rmsd[4], 'b--', label=rdf_sec)
        ax.legend(loc="best")
        ax.set_xlabel("time (ps)")
        ax.set_ylabel(r"RMSD ($\AA$)")
        fig.savefig("rmsd_md_analysis.png")


    def Energychek(self):
        from ase.io import read, write
        import os
        import MDAnalysis
        import MDAnalysis.analysis.rdf
        import MDAnalysis.analysis.rms
        import matplotlib.pyplot as plt
        import pandas as pd
        import re


        os.chdir(self.dire)

        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'8'
        }



        
        df = pd.read_csv('aimd-1.ener',sep='\s+',skipinitialspace=True,engine='python',)


        a1=df['Step']
        a2=df['Nr.']
        a3=df['Time[fs]']
        a4=df['Kin.[a.u.]']
        a5=df['Temp[K]']
        a6=df['Pot.[a.u.]']
        # a7=df['Time[fs]']

        fig = plt.figure(figsize=(3.3464566929*2.1,2.5098425197*2.1),dpi=1200)
        plt.subplot(2,2,1)
        plt.scatter(a1, a2, alpha=0.85, label='Kinetic energy',linewidth=1,c='NONE',edgecolors='blue', s=5)
        # plt. scatter(litt, liss, alpha=0.85, label='Energy',linewidth=2)
        plt.xlabel(r"Time (fs)",fontsize=8,family='Arial')
        plt.ylabel(r"Energy (a.u).",fontsize=8,family='Arial')
        plt.legend(edgecolor='none', prop=font1,)
        plt.tick_params(labelsize=8)

        plt.subplot(2,2,2)
        plt.scatter(a1, a4, alpha=0.85, label='Potential energy',linewidth=1,c='NONE',edgecolors='red', s=5)
        # plt. scatter(litt, liss, alpha=0.85, label='Energy',linewidth=2)
        plt.xlabel(r"Time (fs)",fontsize=8,family='Arial')
        plt.ylabel(r"Energy (a.u).",fontsize=8,family='Arial')
        plt.legend(edgecolor='none', prop=font1,)
        plt.tick_params(labelsize=8)


        plt.subplot(2,2,3) 
        plt.scatter(a1, a5, alpha=0.85, label='Conserved quantity',linewidth=1,c='NONE',edgecolors='grey', s=5)
        # plt. scatter(litt, liss, alpha=0.85, label='Energy',linewidth=2)
        plt.xlabel(r"Time (fs)",fontsize=8,family='Arial')
        plt.ylabel(r"Energy (a.u).",fontsize=8,family='Arial')
        plt.legend(edgecolor='none', prop=font1,)
        plt.tick_params(labelsize=8)

        plt.subplot(2,2,4)
        plt.scatter(a1, a3, alpha=0.85, label='Temperature',linewidth=1,c='NONE',edgecolors='forestgreen', s=5)
        plt.xlabel(r"Time (fs)",fontsize=8,family='Arial')
        plt.ylabel(r"Temp (K)",fontsize=8,family='Arial')
        # plt.ylim((0,10000))
        plt.legend(edgecolor='none', prop=font1,)
        plt.tick_params(labelsize=8)

        # plt.subplot(2,3,5)
        # plt.scatter(a1, a6, alpha=0.85, label='UsedTime',linewidth=2,c='purple')
        # plt.xlabel(r"Time (fs)",fontsize=8,family='Arial')
        # plt.ylabel(r"UsedTime (s)",fontsize=8,family='Arial')
        # plt.legend(edgecolor='none', prop=font1,)
        # plt.tick_params(labelsize=8)
        

        # a8=[]
        # for i, element in enumerate(a6):
        #     nv=a6[:i]
        #     su = sum(nv)
        #     a8.append(su/60)
        
        # plt.subplot(2,3,6)
        # plt.scatter(a1, a8, alpha=0.85, label='Accumtime',linewidth=2,c='cadetblue')
        # # plt.plot(a1, a7, alpha=0.85, label='Energy',linewidth=2)
        # plt.xlabel(r"Time (fs)",fontsize=8,family='Arial')
        # plt.ylabel(r"Accumtime (min)",fontsize=8,family='Arial')
        # plt.legend(edgecolor='none', prop=font1,)
        # plt.tick_params(labelsize=8)



        dm = pd.DataFrame({'Time (fs)':a1,'Kin':a2,'Temp':a3,'Pot':a4,'Cons Qty':a5,'UsedTime':a6})
        dm.to_excel("Energychek.xlsx")


        # 以下是图片的格式设置
        # 设置横纵坐标的名称以及对应字体格式
        font2 = {'family' : 'Arial',
        'weight' : 'normal','size':'x-large'
        }

        # plt.xlim((0,12))
        # plt.ylim((0))
        #不显示Y轴的刻度
        # plt.yticks([])

        #设置图例对应格式和字体
        # auto_adjust_subplotpars(fig, renderer,(2,3) , [1,1,1,1,1,1], (2,3,1), ax_bbox_list=None, pad=1.08, h_pad=None, w_pad=None, rect=None)
        plt.subplots_adjust(wspace =0.3,hspace =0.3 )
        #存储为
        plt.savefig('Energychek.tiff', bbox_inches='tight',transparent=True,format='tiff')#指定分辨率,边界紧，背景透明
        plt.show()
        print ('congratulate！！！')

    def Wullff(self,surface_energy,directions=(3,4,4),x=0.47,y=-16):
        from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice
        # Import the neccesary tools for making a Wulff shape
        from pymatgen.analysis.wulff import WulffShape
        from pymatgen.ext.matproj import MPRester
        from pymatgen.io.cif import CifParser
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        import pandas as pd
        import os
        os.chdir(self.dire)
######################################################################################################################        
        pic_name='wulff'
        structure = CifParser('aimd.cif')
        struct = structure.get_structures()[0]
        structureslab = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
        min_slab_size_1=8.0
        min_vacuum_size_1=15
        
        need_miller_index =surface_energy.keys()#[[i for i in j] for j in surface_energy.keys()]
        Energy_s=surface_energy.values()
        Energy_s=[(float(i)/96.485354)*1.602176462E-19 for i in Energy_s]
        surene=[]
        for i,u in enumerate(need_miller_index):         
            slab = SlabGenerator(struct, miller_index=u, min_slab_size=min_slab_size_1,\
                                 min_vacuum_size=min_vacuum_size_1, center_slab=True)        
            slabs_bak = slab.get_slabs()[0].copy()#可能的晶面
            # slabs.make_supercell(self.supercell)
            #晶胞扩充
            cc = slabs_bak.surface_area*1E-20
            print ('The surface area of the plane %s is %s square meter\n\n'%(u,cc))

            print ('The surface tension of the plane %s is %s Joule per square meter \n\n'%(u,Energy_s[i]/cc))
            surene.append(Energy_s[i]/cc)

        for i,j in enumerate(need_miller_index):
            surface_energy[j]=surene[i]

        print('###################################################################')
        print(surface_energy)
        print('###################################################################')

######################################################################################################################        
# wulff
        surface_energies = surface_energy
        miller_list = surface_energies.keys()
        e_surf_list = surface_energies.values()
        font2 = {'family' : 'Times New Roman','fontsize':'40' , 'weight' : 'bold',}
        wulffshape = WulffShape(struct.lattice, miller_list, e_surf_list)
        print(wulffshape.area_fraction_dict)
        os.chdir(self.dire)
        dict1 = wulffshape.area_fraction_dict
        xx =[]
        yy =[]
        for key,value in dict1.items():
            xx.append(key)
            yy.append(value)
 
        cc=wulffshape.effective_radius
        bb =wulffshape.volume
        dd = wulffshape.shape_factor
        print(wulffshape.effective_radius)
        res=pd.DataFrame({'Slab':xx,'area':yy})#构造原始数据文件
        df = res.sort_values(by='area', ascending=True) 
        with open(str(pic_name)+'.txt','w') as f:

            f.write('effective radius:' +str(cc)+'         '+"volume:"+str(bb)+'         '+"shape factor:"+str(dd))
        print (cc)       
        # os.chdir(r"D:\Desktop\VASP practical\workdir")
        df.to_excel(str(pic_name)+".xlsx")#生成Excel文件，并存到指定文件路径下
        wulffshape.get_plot(bar_on=True,aspect_ratio=(8,8) ,bar_pos=[0, 0.85, 1.1, 0.045],direction=directions)


        plt.title(str(pic_name),font2,x=x,y=y)  
        
        plt.savefig(str(pic_name)+".tiff",bbox_inches='tight', transparent=True,dpi=600,format='tiff')
        
    def Elec_poten_peri(self,filename='aimd-pos-1.pdb'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import fnmatch
        import subprocess
        os.chdir(self.dire)

######################################################################################################################
        filename01 = str('aimd-v_hartree-1_0.cube')
######################################################################################################################
        inputli=['13','0','ESP1.cub','-1','q']
        filenas01 = str('cubchange01.txt')
        sg='\n'.join(inputli)
        fp=open(filenas01,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
        filename02 = str('aimd-ELECTRON_DENSITY-1_0.cube')
######################################################################################################################
        inputli=['13','0','density1.cub','-1','q']
        filenas02 = str('cubchange02.txt')
        sg='\n'.join(inputli)
        fp=open(filenas02,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd_01='Multiwfn %s < cubchange01.txt\n'%str(filename01)
        cmd_02='Multiwfn %s < cubchange02.txt\n'%str(filename02)
        cmd_04='copy /Y *.cub %s\n'%str("D:\Software\VMD")
        res=os.popen(cmd_01)
        output_str=res.read() 
        print(output_str) 
        res=os.popen(cmd_02)
        output_str=res.read() 
        print(output_str)   
        os.remove(filename01)
        os.remove(filename02)
        res=os.popen(cmd_04)
        output_str=res.read()
        print(output_str)              
        os.remove(filenas01)
        os.remove(filenas02)






    def Elec_potential(self,filename='aimd-pos-1.pdb'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import fnmatch
        import subprocess
        os.chdir(self.dire)

######################################################################################################################
        filename = str('molden_for_multiwfn.molden')
######################################################################################################################
        inputli=['5','1','2','2','0','5','12','1','2','0','r','totesp.cub','13','11','5','27.2114','0','totesp.cub']
        filenas = str('ESPiso.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd_01='Multiwfn %s < ESPiso.txt\n'%str(filename)
        cmd_02='move /Y density.cub density1.cub\n'
        cmd_03='move /Y totesp.cub ESP1.cub\n'
        cmd_04='move /Y *.cub %s\n'%str("D:\Software\VMD")
        res=os.popen(cmd_01)
        output_str=res.read() 
        print(output_str) 
        res=os.popen(cmd_02)
        output_str=res.read() 
        print(output_str)   
        res=os.popen(cmd_03)
        output_str=res.read()
        print(output_str)
        res=os.popen(cmd_04)
        output_str=res.read()
        print(output_str)              
        os.remove(filenas)

######################################################################################################################
        inputli=['12','3','0.15','0','5','mol.pdb','6']
        filenas = str('ESPiso.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd_01='Multiwfn %s < ESPpt.txt\n'%str(filename)
        cmd_02='move /Y vtx.pdb vtx1.pdb\n'
        cmd_03='move /Y mol.pdb mol1.pdb\n'
        cmd_04='copy /Y *.pdb %s\n'%str("D:\Software\VMD")
        res=os.popen(cmd_01)
        output_str=res.read()
        print(output_str)
        res=os.popen(cmd_02)
        output_str=res.read() 
        print(output_str)   
        res=os.popen(cmd_03)
        output_str=res.read()
        print(output_str)
        res=os.popen(cmd_04)
        output_str=res.read()
        print(output_str)              
        os.remove(filenas)
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')

######################################################################################################################
        inputli=['12','3','0.15','0','5','mol.pdb','6','2']
        filenas = str('ESPext.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd_01='Multiwfn %s < ESPext.txt\n'%str(filename)
        cmd_02='move /Y vtx.pdb vtx1.pdb\n'
        cmd_03='move /Y mol.pdb mol1.pdb\n'
        cmd_04='copy /Y *.pdb %s\n'%str("D:\Software\VMD")
        res=os.popen(cmd_01)
        output_str=res.read()
        print(output_str)
        res=os.popen(cmd_02)
        output_str=res.read() 
        print(output_str)   
        res=os.popen(cmd_03)
        output_str=res.read()
        print(output_str)
        res=os.popen(cmd_04)
        output_str=res.read()
        print(output_str)              
        os.remove(filenas)
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')
######################################################################################################################

    def vmdrender(self,):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import fnmatch
        import subprocess
        
        os.chdir(self.dire)
######################################################################################################################
######################################################################################################################
        cmd_01=r"copy /Y D:\Software\VMD\vmdscene.dat vmdscene.dat"
        cmd_02='tachyon_WIN32 vmdscene.dat -format BMP -o VMD.bmp -aasamples 24 -mediumshade -trans_vmd -res 2000 1500 -numthreads 4\n'
        res=os.popen(cmd_01)
        output_str=res.read() 
        print(output_str)              
        res=os.popen(cmd_02)
        output_str=res.read() 
        print(output_str) 
      
######################################################################################################################



    def pdbtocif(self,filename='aimd-pos-1.pdb'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import fnmatch
        os.chdir(self.dire)
######################################################################################################################
        for root,dirs,files in os.walk(os.getcwd()):
            if filename in files:
                print(root)
                # print(files)
                filename='aimd-pos-1.pdb'
            else:
                filename='aimd.pdb'
######################################################################################################################
        f=open(filename,'r')
        alllines=f.readlines()#.strip().split()
        f.close()
        indexs=[n for n, u in enumerate(alllines) if 'CRYST1' in u ]
        alllines=alllines[:1]+alllines[indexs[-1]:]
######################################################################################################################
#WRITE
        filename = str('aimd.pdb')
        s=''.join(alllines)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()
        print(fp.close())
        print('The last frame is extracted')        
######################################################################################################################
        inputli=['\n','100','2','33','\n','0','q']
        filenas = str('tes.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd='Multiwfn %s -isilent < tes.txt'%str(filename)
        # print(cmd)
        res=os.popen(cmd)
        output_str=res.read()
        # print(output_str)
        os.remove(filenas)
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')
######################################################################################################################
        for root,dirs,files in os.walk(os.getcwd()):
            if 'aimd-pos-1.pdb' in files:
                print(files)
                os.remove('aimd-pos-1.pdb')
                print('    You have rewritten aimd-pos-1.pdb to aimd.cif            ')
            else:
                print('          You have rewritten aimd.pdb to aimd.cif        ')
          



    def forcechek(self,atomindex):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        os.chdir(self.dire)
        ######################################################################################################################
        '''
        提取数据
        '''
        filename='aimd-frc-1.xyz'
        def redcsv(filenamelist=filename):#read file for eachline
            f=open(filenamelist,'r')
            alllines=f.readlines()
            f.close()
            Atomsnumber=[]
            Step=[]
            Energy=[]
            Forcesx=[]
            Forcesy=[]
            Forcesz=[]

            for i,eachline in enumerate(alllines):
                if i == 0:
                    Atomsnumber.append(float(eachline.split()[-1]))
            for i,eachline in enumerate(alllines):
                if 'E =' in eachline:
                    Step.append(float(eachline.split()[2][:-1]))
                    Energy.append(float(eachline.split()[-1]))

                    inum=atomindex # please input the number of atom

                    Forcesx.append(float(alllines[i+inum].split()[-3]))
                    Forcesy.append(float(alllines[i+inum].split()[-2]))
                    Forcesz.append(float(alllines[i+inum].split()[-1]))


            return [Atomsnumber, Step,Energy,Forcesx,Forcesy,Forcesz]
        ######################################################################################################################
        '''
        绘图
        '''
        dvale=redcsv()
        xss=dvale[1]#step
        yfe=dvale[2]#energy
        yfe=[round(i,1) for i in yfe]
        yfx=dvale[3]#x_force
        yfy=dvale[4]#y_force
        yfz=dvale[5]#z_force

        fig,(ax1, ax2) = plt.subplots(1,2,figsize=(17.5*self.inch_cm,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)
        ax1.plot(xss,yfe,alpha=0.75,linewidth=2.8,label='Energy')
        ax2.plot(xss,yfx,alpha=0.75,linewidth=2.4,label='Force component X')
        ax2.plot(xss,yfy,alpha=0.75,linewidth=2.4,label='Force component Y')
        ax2.plot(xss,yfz,alpha=0.75,linewidth=2.4,label='Force component Z')
        # plt.axhline(y=2,c="black")
        # plt.axvline(x=534,c="green")#461-500,520-534,
        # 以下是图片的格式设置
        # 设置横纵坐标的名称以及对应字体格式
        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'8'
        }
        # ax.legend(prop=font1,frameon=False,)
        # ax.set_xlim((320,600))
        ax1.set_xlabel(r'Step',fontsize=8,family='Arial')
        ax1.set_ylabel(r'Energy (a.u.)',fontsize=8,family='Arial')
        ax2.set_xlabel(r'Step',fontsize=8,family='Arial')
        ax2.set_ylabel(r'Force (Hartree·bohr$^{-1}$)',fontsize=8,family='Arial')

        # ax1 = ax.twiny()
        # ax1.scatter(times, weight02,marker='o',alpha=0.95,color='NONE', s=5,edgecolors='blue',linewidth=0.8,label='Weight-Time')
        # font1 = {'family' : 'Arial',
        # 'weight' : 'normal','size':'8'
        # }
        # ax1.legend(prop=font1,frameon=False,)
        # ax.set_xlim((-10,270))
        # ax1.set_xlabel(r'Time (min)',fontsize=8,family='Arial')
        # ax1.set_ylabel(r'Weight (%)',fontsize=8,family='Arial')
        # cd=max(a5)/10
        # ax.set_xlim((max(wavenumber),min(wavenumber)))
        # ax1.set_xlim((20,40))
        #不显示Y轴的刻度
        plt.legend(prop=font1,frameon=False,ncol=1)
        # plt.xticks([i for i in np.arange(0,261,52)])
        # plt.yticks([i for i in np.arange(0,1.01,0.25)])
        # ax1.set_yticks([])
        # ax2.set_yticks([])
        #存储为
        os.chdir(self.dire)
        plt.savefig('forcechek.tiff', transparent=True,format='tiff')#指定分辨率,边界紧，背景透明bbox_inches='tight',
        plt.show()
        print ('congratulate！！！')


    def freeenergy(self):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        os.chdir(self.dire)
        ######################################################################################################################
        '''
        提取数据
        '''
        bohr_2_angstrom = 0.529177
        kb = 8.6173303e-5 # eV * K^-1

        temperature = 1000.0                          # Change the temperature according to your MD simulations!
        colvar_path = "aimd-COLVAR.metadynLog"

        # Load the colvar file
        colvar_raw = np.loadtxt(colvar_path)

        # Extract the two CVs
        d1 = colvar_raw[:, 1] * bohr_2_angstrom
        d1min=min(d1)
        d1max=max(d1)
        d2 = colvar_raw[:, 2] * bohr_2_angstrom
        d2min=min(d2)
        d2max=max(d2)
        # print(colvar_raw[:, 2])

        # Create a 2d histogram corresponding to the CV occurances
        binvalue=[120,120]
        cv_hist = np.histogram2d(d1, d2, bins=binvalue)
        # print(cv_hist[0])

        # probability from the histogram
        prob = cv_hist[0]/len(d1)
        # print(prob)

        # Free energy surface
        fes = -kb * temperature * np.log(prob)
        chek_where=np.isinf(fes)
        fes[chek_where]=0
        ######################################################################################################################
        '''
        绘图
        '''
        xmin=d1min
        ymin=d2min
        xmax=d1max
        ymax=d2max
        intensity=fes
        # print(intensity.shape)
        g_conv, δ_conv = np.mgrid[xmin:xmax:(xmax-xmin)/(binvalue[0]),ymin:ymax:(ymax-ymin)/(binvalue[1])]


        maxint=max([max(i) for i in intensity])
        # print(maxint)
        intensity[intensity==0]=maxint
        minint=min([min(i) for i in intensity])
        cbarticks=np.arange(round(minint, 1),round(maxint, 1),(round(maxint, 1)-round(minint, 1))/5)
        # print((round(minint, 2)-round(maxint, 2))/5)
        # levels = np.arange(0, maxint, maxint*0.1)
        # points=np.meshgrid(wavenumber,time)
        # points=[list(list(i) for i in i) for i in points]
        # print(intensity)


        fig,ax = plt.subplots(figsize=(17.5*self.inch_cm*0.55,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)
        


        cse = ax.contourf(g_conv, δ_conv,intensity,alpha=0.65,cmap ='terrain',levels=200)#cmap ="gist_heat",binary
        # contour = ax.contour(x1, y1,deltat,[0],colors=c,linestyles='-',label='d')
         
        # contour = ax.contour(g_conv, δ_conv,intensity,levels,linewidths=0.5, colors='k')
        # plt.axhline(y=1, xmin=0.0, c='grey',ls='--',lw=1) 
        # plt.axvline(x=24.92, ymin=0,  c='red',ls='-.',lw=0.8) 
        # plt.axvline(x=14.06, ymin=0,  c='red',ls='-.',lw=0.8)
        # plt.axvline(x=29.15, ymin=0,  c='red',ls='-.',lw=0.8)  
        ax.set_xlabel(r"R$_i$ ($\.A$)",fontsize=8,family='Arial')
        ax.set_ylabel(r"R$_j$ ($\.A$)",fontsize=8,family='Arial')
        # ax.legend(edgecolor='none',facecolor='none', prop=font1,loc='best')
        ax.tick_params(labelsize=8)
        cbar = fig.colorbar(cse)
        cbar.ax.set_ylabel("Free energy (eV)",fontsize=8,family='Arial')
        cbar.set_ticks(cbarticks)
        cbar.ax.tick_params(labelsize=8)
        # Add the contour line levels to the colorbar
        # cbar.add_lines(cse)
        # plt.colorbar(cse)
        # 以下是图片的格式设置
        # 设置横纵坐标的名称以及对应字体格式
        font2 = {'family' : 'Arial',
        'weight' : 'normal','size':'x-large'
        }

        # ax.set_ylim((0,250))
        # cd=max(a5)/10
        # ax.set_xlim((max(wavenumber),min(wavenumber)))
        # ax.set_xlim((10,50))
        #不显示Y轴的刻度
        # plt.xticks([i for i in np.arange(-0.2,-0.099,0.05)])
        # plt.yticks([i for i in np.arange(-0.2,0.21,0.1)])
        # ax.set_yticks([])
        #存储为
        os.chdir(self.dire)
        plt.savefig('./part%s.tiff'%str(self.dire[-2:]),transparent=True,format='tiff')#指定分辨率,边界紧，背景透明 bbox_inches='tight'
        plt.show()
        print ('##############         DOWN  congratulate！！！  #######################')



    def pdos(self,nelement):
        import os
        import matplotlib.pyplot as plt
        import numpy as np
        os.chdir(self.dire)
        ######################################################################################################################
        '''
        提取数据
        '''
        # Parameters
        FWHM = 0.4   # Full width half maximum
        Emin=-16.    # Starting energy
        Emax= 4.     # -----
        DE = 0.01    # Fineness of the energies to convolve to
        dE=np.arange(Emin,Emax,DE)
        E_shift = 0
        filename = str('aimd_pdos-k%s-1.pdos'%str(nelement))

        input_file = open(filename, 'r')

        firstline  = input_file.readline().strip().split()
        secondline = input_file.readline().strip().split()

        # Kind of atom
        atom = firstline[6]
        print(atom)
        #iterationstep
        iterstep = int(firstline[12][:-1]) #[:-1] delete ","
        # Energy of the Fermi level
        efermi = float(firstline[15])
        print(efermi)

        if 'd' in secondline:            
            def convolve(E,PDOS):
                
                sigma = 1 / (2 * np.sqrt(2 * np.log(2) ) ) * FWHM

                conv = np.zeros((len(dE),3))
                for j in range(conv.shape[1]):
                    W = PDOS[:,j]
                    for i,e in enumerate(E):
                        conv[:,j] += np.exp(-(dE-e)**2 / (2 * sigma**2))*W[i]
                return conv[:,0], conv[:,1],conv[:,2]
        else:
            def convolve(E,PDOS):
                
                sigma = 1 / (2 * np.sqrt(2 * np.log(2) ) ) * FWHM

                conv = np.zeros((len(dE),2))
                for j in range(conv.shape[1]):
                    W = PDOS[:,j]
                    for i,e in enumerate(E):
                        conv[:,j] += np.exp(-(dE-e)**2 / (2 * sigma**2))*W[i]
                return conv[:,0], conv[:,1]            



                
        PDOS = np.asarray(open(filename).read().split())

        if 'f' in PDOS:
            start = 26
            Norb = 4+3
        elif 'd' in PDOS:
            start = 25
            Norb = 3+3
        elif 'p' in PDOS:
            start = 24
            Norb = 2+3
        else:
            start = 23
            Norb = 1+3
            
        L = int(PDOS[-Norb])
        
        pdos = PDOS[start:]
        # print(pdos)
        pdos = np.asarray(np.split(pdos,L))
        # print (pdos)
        E = pdos[:,1].astype(np.float)*27.2114 #a.u. -> eV
        # print (pdos[:,3:])
        pdos = pdos[:,3:].astype(float)

        if 'd' in secondline:
            s, p,d= convolve(E,pdos)
            tdos=np.array(s) + np.array(p)+ np.array(d)
        else:
            s, p= convolve(E,pdos)
            tdos= np.array(s) + np.array(p)


        ######################################################################################################################
        '''
        绘图
        '''
        fig,ax = plt.subplots(figsize=(17.5*self.inch_cm*0.55,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)
        ax.fill_between(dE+E_shift,tdos, alpha=.75,lw=.25,zorder=2,facecolor='#CCCCCC',label ='TDOS')
        ax.plot(dE+E_shift,s, c='blue',linestyle='-',linewidth=1., label='s')
        # ax.fill_between(dE+E_shift,s, alpha=.5,lw=.25,zorder=2,facecolor='#663333',label='%s s'%str(atom))
        ax.plot(dE+E_shift,p, c='red',linestyle='-',linewidth=1., label = 'p')
        # ax.fill_between(dE+E_shift,p, alpha=.5,lw=.25,zorder=2,facecolor='#339999',label = '%s p'%str(atom))
        
        if 'd' in secondline:
            ax.plot(dE+E_shift,d, c='black', alpha=1,linestyle='-',linewidth=1., label = 'd')
            # ax.fill_between(dE+E_shift,d, alpha=.5,lw=.25,zorder=2,facecolor='#CCCC66',label = '%s d'%str(atom))
        # ax.plot(dE+E_shift,tdos, c='blue',linestyle='-',linewidth=1., label = '%s TDOS'%str(atom))
        

        plt.axvline(x=efermi,  c='grey',ls='dashdot',lw=1)


        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'8'
        }
        # ax.set_ylim((-3,80))
        # ax.set_xlim((0,250))
        ax.legend(prop=font1,frameon=False,ncol=1)
        ax.tick_params(labelsize=8)        
        ax.set_xlabel('Binding energy (eV)',fontsize=8,family='Arial')
        ax.set_ylabel('PDOS',fontsize=8,family='Arial')
        # ax.legend(edgecolor='none',facecolor='none', prop=font1,loc='best')

        ax.set_yticks([])
        # ax.set_xticks([])
        ax.tick_params(labelsize=8)
        # plt.yticks([])                                                                                              
        # plt.ylabel('PDOS (arbitrary units)')
        # plt.xlabel('Binding energy (eV)')
        # plt.xticks(np.arange(-15,5,5),-np.arange(-15,5,5))
        # plt.xlim([-16,3])
        # plt.title('FWHM = ' + str(FWHM) + ' eV')

        # plt.legend(loc = 1)

        plt.savefig('pdos.tiff', bbox_inches='tight',transparent=True,format='tiff')#指定分辨率,边界紧，背景透明
        plt.show()
        print ('congratulate！！！')

    def lpdos(self,nlist):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        os.chdir(self.dire)
        ######################################################################################################################
        '''
        提取数据
        '''
        # Parameters
        FWHM = 0.5   # Full width half maximum
        Emin=-16.    # Starting energy
        Emax= 4.     # -----
        DE = 0.01    # Fineness of the energies to convolve to
        dE=np.arange(Emin,Emax,DE)
        E_shift = 0
        filename = str('aimd_pdos-list%s-1.pdos'%str(nlist))

        input_file = open(filename, 'r')
        
        firstline  = input_file.readline().strip().split()
        # print(firstline)
        secondline = input_file.readline().strip().split()
        # print(secondline)

        #Kind of atom
        # atom = firstline[6]
        # print(atom)
        #iterationstep
        # iterstep = int(firstline[12][:-1]) #[:-1] delete ","
        # Energy of the Fermi level
        efermi = float(firstline[-2])
        print(efermi)

        if 'd' in secondline:            
            def convolve(E,PDOS):
                
                sigma = 1 / (2 * np.sqrt(2 * np.log(2) ) ) * FWHM

                conv = np.zeros((len(dE),3))
                for j in range(conv.shape[1]):
                    W = PDOS[:,j]
                    for i,e in enumerate(E):
                        conv[:,j] += np.exp(-(dE-e)**2 / (2 * sigma**2))*W[i]
                return conv[:,0], conv[:,1],conv[:,2]
        else:
            def convolve(E,PDOS):
                
                sigma = 1 / (2 * np.sqrt(2 * np.log(2) ) ) * FWHM

                conv = np.zeros((len(dE),2))
                for j in range(conv.shape[1]):
                    W = PDOS[:,j]
                    for i,e in enumerate(E):
                        conv[:,j] += np.exp(-(dE-e)**2 / (2 * sigma**2))*W[i]
                return conv[:,0], conv[:,1]            



                
        PDOS = np.asarray(open(filename).read().split())
        print(PDOS)

        if 'f' in PDOS:
            start = 28
            Norb = 4+3
        elif 'd' in PDOS:
            start = 27
            Norb = 3+3
        elif 'p' in PDOS:
            start = 26
            Norb = 2+3
        else:
            start = 25
            Norb = 1+3
            
        L = int(PDOS[-Norb])
        
        pdos = PDOS[start:]
        print(pdos)
        pdos = np.asarray(np.split(pdos,L))
        # print (pdos)
        E = pdos[:,1].astype(np.float)*27.2114 #a.u. -> eV
        # print (pdos[:,3:])
        pdos = pdos[:,3:].astype(float)
        if 'd' in secondline:
            s, p,d= convolve(E,pdos)
            tdos=np.array(s) + np.array(p)+ np.array(d)
        else:
            s, p= convolve(E,pdos)
            tdos= np.array(s) + np.array(p)


        ######################################################################################################################
        '''
        绘图
        '''
        fig,ax = plt.subplots(figsize=(17.5*self.inch_cm*0.55,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)
        ax.plot(dE+E_shift,s, c='blue',linestyle='-',linewidth=1., label='s')
        # ax.fill_between(dE+E_shift,s, alpha=.5,lw=.25,zorder=2,facecolor='#663333',label='%s s'%str(atom))
        ax.plot(dE+E_shift,p, c='red',linestyle='-',linewidth=1., label = 'p')
        # ax.fill_between(dE+E_shift,p, alpha=.5,lw=.25,zorder=2,facecolor='#339999',label = '%s p'%str(atom))
        if 'd' in secondline:
            ax.plot(dE+E_shift,d, c='black',linestyle='-',linewidth=1., label = 'd')
            # ax.fill_between(dE+E_shift,d, alpha=.5,lw=.25,zorder=2,facecolor='#CCCC66',label = '%s d'%str(atom))
        # ax.plot(dE+E_shift,tdos, c='blue',linestyle='-',linewidth=1., label = '%s TDOS'%str(atom))
        ax.fill_between(dE+E_shift,tdos, alpha=.75,lw=.25,zorder=2,facecolor='#CCCCCC',label ='TDOS')
        plt.axvline(x=efermi,  c='grey',ls='dashdot',lw=1)


        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'8'
        }
        # ax.set_ylim((-3,80))
        # ax.set_xlim((0,250))
        ax.legend(prop=font1,frameon=False,ncol=1)
        ax.tick_params(labelsize=8)        
        ax.set_xlabel('Binding energy (eV)',fontsize=8,family='Arial')
        ax.set_ylabel('PDOS',fontsize=8,family='Arial')
        # ax.legend(edgecolor='none',facecolor='none', prop=font1,loc='best')

        ax.set_yticks([])
        # ax.set_xticks([])
        ax.tick_params(labelsize=8)
        # plt.yticks([])                                                                                              
        # plt.ylabel('PDOS (arbitrary units)')
        # plt.xlabel('Binding energy (eV)')
        # plt.xticks(np.arange(-15,5,5),-np.arange(-15,5,5))
        # plt.xlim([-16,3])
        # plt.title('FWHM = ' + str(FWHM) + ' eV')

        # plt.legend(loc = 1)

        plt.savefig('pdos.tiff', bbox_inches='tight',transparent=True,format='tiff')#指定分辨率,边界紧，背景透明
        plt.show()
        print ('congratulate！！！')     

    def IRspectrum(self,corr=100,FWHM=15,xlim_1=800,xlim_2=1200):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        os.chdir(self.dire)
        def minmaxx(aa):#归一化到最大值

            lisnew=[(i-min(aa))/(max(aa)-min(aa)) for i in aa]

            return lisnew
        ######################################################################################################################
        '''
        提取数据
        '''
        # Parameters
        FWHM = FWHM # Full width half maximum
        Emin=300.    # Starting energy
        Emax=4000.     # -----
        DE = 1    # Fineness of the energies to convolve to
        dE=np.arange(Emin,Emax,DE)
        # dE=minmaxx(dE)

        def convolve(FREQ,INT):


            # sigma = 1 / (2 * np.sqrt(2 * np.log(2) ) ) * FWHM

            conv = np.zeros((len(dE),1))

            for n,j in enumerate(dE):
                for u,v in enumerate(FREQ):
                    if dE[n]<v<dE[n+1]:
                        conv[n]=INT[u]

            cont=[]
            nozero=[i for i,j in enumerate(conv) if j !=0 ]
            
            for u,v in enumerate(conv):
                hus=[]
                for i in nozero:
                    hu=(np.exp(-(u-dE[i])**2 / (2 *FWHM**2))*conv[i])
                    hus.append(hu)
                huus=sum(hus)

                cont.append(huus)
            return dE,cont 


        filename = str('aimd-VIBRATIONS-1.mol')
        f=open(filename,'r')
        alllines=f.readlines()
        f.close()
        for i,eachline in enumerate(alllines):
            if '[FREQ]' in eachline:
                state=i+1
            elif '[FR-COORD]' in eachline:
                end=i
            elif '[INT]' in eachline:
                began=i+1
        cellvetor_1=alllines[state:end]
        cellvetor_2=alllines[began:]
        cellvetor_1=[i for i in cellvetor_1 if not "*" in i]
        cellvetor_2=[cellvetor_2[n] for n,i in enumerate(cellvetor_1) if not "*" in i ]

        cellvetor_1=[float(i) for i in cellvetor_1 if 4000>float(i)>0]
        cellvetor_2=[float(cellvetor_2[n]) for n,i in enumerate(cellvetor_1) if 4000>float(i)>0]

        # print (cellvetor_1)
        FREQ= cellvetor_1
        FREQ=[i-corr for i in FREQ]
        INT=cellvetor_2
        INT=minmaxx(INT)



        freq = np.asarray(FREQ)
        intensity = np.asarray(INT)

        intensity= convolve(freq,intensity)

        # # print(intensity[0],intensity[1])
        # ######################################################################################################################
        # '''
        # 绘图
        # '''
        fig,ax = plt.subplots(figsize=(17.5*self.inch_cm,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)
        ax.plot(intensity[0],intensity[1], c='blue',linestyle='-',linewidth=1., label='Calculated')
        plt.axvline(x=xlim_1,  c='red',ls='dashdot',lw=1)
        plt.axvline(x=xlim_2,  c='purple',ls='dashdot',lw=1)

        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'8'
        }
        # ax.set_ylim((-3,80))
        ax.set_xlim((4000,300))
        ax.legend(prop=font1,frameon=False,ncol=1)
        ax.tick_params(labelsize=8)        
        ax.set_xlabel(r'Wavenumber (cm $^{-1}$)',fontsize=8,family='Arial')
        ax.set_ylabel('Intensity (a.u.)',fontsize=8,family='Arial')
        # ax.legend(edgecolor='none',facecolor='none', prop=font1,loc='best')

        ax.set_yticks([])
        # ax.set_xticks([])
        ax.tick_params(labelsize=8)

        plt.savefig('irspectrum.tiff', bbox_inches='tight',transparent=True,format='tiff')#指定分辨率,边界紧，背景透明
        plt.show()
        print ('congratulate！！！')     

    def moldenformultiwfn(self):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        import re
        os.chdir(self.dire)
######################################################################################################################
        os.chdir(r'F:\Computing tutorial\cp2k\cp2k-7.1\cp2k-7.1\data')
        filename = str('POTENTIAL')
        f=open(filename,'r')
        alllines=f.readlines()
        f.close()  
        os.chdir(self.dire)   

        elemlist=[]   
        Valence_electron_number=[]
        for i,eachline in enumerate(alllines):
            if 'GTH-PBE-q' in eachline: 
                elel=eachline.strip().split()[0]
                elemlist.append(elel)
                val=eachline.strip().split()[1][9:]
                Valence_electron_number.append(float(val))
        # print(Valence_electron_number)


        ven=Valence_electron_number


######################################################################################################################
        #Read cell vector
        filename = str('aimd.inp')
        f=open(filename,'r')
        alllines=f.readlines()
        f.close()
        for i,eachline in enumerate(alllines):
            if '&CELL' in eachline:
                state=i
        cellvetor=alllines[state+1:state+4]
        cellvetor=[i.strip().split() for i in cellvetor]
        cellvetor=np.array(cellvetor).reshape((3,4))
        cellvetor=cellvetor[:,1:]
        natomsv=open('aimd.pdb').readline().strip().split()
        natoms=int(natomsv[-2])

        element=[]
        for i in alllines:
            if "ELEMENT" in i:
                eles=i.strip().split()[-1]
                element.append(eles)
        element=list((set(element)))
        print('all the elements are %s'%str(element))
        # print('The namber of atoms is %s'%str(natoms))


# ######################################################################################################################
        filename = str('aimd-MOS-1_0.molden')
        fp=open(filename)
        lines=[]
        for line in fp:
            lines.append(line)
        fp.close()
        #insert lattice vector
        lines.insert(1,' [Cell]'+'\n')
        lines.insert(2,' '+'     '.join(cellvetor[0])+'\n')
        lines.insert(3,' '+'     '.join(cellvetor[1])+'\n')
        lines.insert(4,' '+'     '.join(cellvetor[2])+'\n')

        #insert valence electron number
        lines.insert(5,' [Nval]'+'\n')
        for n,i in enumerate(element):
            lines.insert(n+6,' %s %s'%(i,ven[elemlist.index(i)])+'\n')
        # lines.insert(7,' H 1'+'\n')
        # lines.insert(8,' S 6'+'\n')        
        # lines.insert(9,' O 6'+'\n')
        # lines.insert(9,' O 6'+'\n')

        s=''.join(lines)
        fp=open('molden_for_multiwfn.molden','w+')
        fp.write(s)
        fp.close()
        print(' Please find your wave function file in the current directory ')




    def IGA(self,):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import fnmatch
        os.chdir(self.dire)
######################################################################################################################
        filename='molden_for_multiwfn.molden'  
######################################################################################################################
        inputli=['\n','20','1','9','\n','0.15\n','3','0','q']
        filenas = str('tes.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
        Number_of_cores=64
        Number_of_nodes=Number_of_cores//64
        inputline=[]
        line1_1='#!/bin/bash'
        inputline.append(line1_1)
        line1_2='#SBATCH -p amd_256 '
        inputline.append(line1_2)  
        line2_1='#SBATCH -N %s'%str(Number_of_nodes)
        inputline.append(line2_1)
        line2_1='#SBATCH -n %s'%str(Number_of_cores)
        inputline.append(line2_1)
        line2_2='#SBATCH -J wcl'
        inputline.append(line2_2)
        line2_2='source /public1/soft/modules/module.sh'
        inputline.append(line2_2)
        line2_2='module load mpi/openmpi/4.1.1'
        inputline.append(line2_2)
        line2_2='export PATH=/public1/soft/openmpi/openmpi-4.1.1-icc/bin:$PATH'
        inputline.append(line2_2)     
        line2_2='export LD_LIBRARY_PATH=/public1/soft/openmpi/openmpi-4.1.1-icc/lib:$LD_LIBRARY_PATH'
        inputline.append(line2_2)
        line2_2='export LIBRARY_PATH=/public1/soft/openmpi/openmpi-4.1.1-icc/lib:$LIBRARY_PATH'
        inputline.append(line2_2)   
        line2_5='Multiwfn %s -isilent < tes.txt'%str(filename)
        inputline.append(line2_5)
        print('                 The number of cores is %s'%str(Number_of_cores))

######################################################################################################################
#WRITE
        filename = str('sub.sh')

        s='\n'.join(inputline)
        fp=open(filename,'wb')
        fp.write(s.encode('utf-8'))
        fp.close()

######################################################################################################################
        print('Please submit the sub file to your supercomputer ')
          

    def localizecube(self,filename='func1.cub',vdw_atom=4,atomindex=83):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import fnmatch
        os.chdir(self.dire)
######################################################################################################################
          
######################################################################################################################

        inputli=['\n','13','13','%s\n'%str(vdw_atom),'20000\n','2','%s\n'%str(atomindex),'0','%s_loc.cub\n'%str(filename[:-4]),'-1','q']
        filenas = str('tes.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
######################################################################################################################
#perform Multiwfn
        cmd='Multiwfn %s -isilent < tes.txt'%str(filename)
        # print(cmd)
        res=os.popen(cmd)
        output_str=res.read()
        # print(output_str)
        os.remove(filenas)
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')
######################################################################################################################
        print('Please find your cub file in current file ')
          

    def den_dif(self,vdw_atom=4,atomindex=83):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import fnmatch
        os.chdir(self.dire)
######################################################################################################################
        filename='minuend.cube'  
######################################################################################################################

        inputli=['\n','13','11','4','subtraction.cube\n','13','%s\n'%str(vdw_atom),'0\n','2','%s\n'%str(atomindex),'0','EDD.cub\n','-1','q']
        filenas = str('tes.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
######################################################################################################################
#perform Multiwfn
        cmd='Multiwfn %s -isilent < tes.txt'%str(filename)
        # print(cmd)
        res=os.popen(cmd)
        output_str=res.read()
        # print(output_str)
        os.remove(filenas)
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')
######################################################################################################################
        print('Please find your cub file in current file ')



    def generateldosinp(self):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        import re
        os.chdir(self.dire)
######################################################################################################################


######################################################################################################################
        #Read cell vector
        filename = str('aimd.inp')
        f=open(filename,'r')
        alllines=f.readlines()
        f.close()
        for i,eachline in enumerate(alllines):
            if '&CELL' in eachline:
                state=i
            elif '&COORD' in eachline:
                began=i+1
            elif '&END COORD' in eachline:
                end=i
        cellvetor=alllines[state+1:state+4]
        cellvetor=[i.strip().split() for i in cellvetor]
        cellvetor=np.array(cellvetor).reshape((3,4))
        cellvetor=cellvetor[:,1:]
        natoms=end-began
        print("############  The number of atoms is %s ###########"%str(natoms))


# ######################################################################################################################


        filename = str('lpdos.inp')


        lines=['&PRINT','    &PDOS','      FILENAME ./aimd_pdos']

        for i in np.arange(natoms):
            strline='     &LDOS\n'+'          LIST %s\n'%str(i)+'     &END LDOS'
            lines.append(strline)

        lines.append('  &END PDOS')
        lines.append('&END PRINT')

        s='\n'.join(lines)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()


    def Ensemblechek(self,filename):
        from ase.io import read, write
        import os
        import MDAnalysis
        import MDAnalysis.analysis.rdf
        import MDAnalysis.analysis.rms
        import matplotlib.pyplot as plt
        import pandas as pd
        import re
        import numpy as np

        os.chdir(self.dire)

        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'24'
        }

        f=open(filename,'r')
        alllines=f.readlines()
        f.close()
        pressurelis=[]
        for i in alllines:
            if 'Pressure' in i:
                pressurelist = i.split()
                pressurelis.append(float(pressurelist[-1]))



        steplis=[]
        steplis.append(0)
        for i in alllines:
            if 'Step number' in i:
                steplist = i.split()
                steplis.append(steplist[-1])
        pressureint=[float(i)*1000//(1000*1.01325) for i in pressurelis ][1:]


        cell_volumelis=[]
        for i in alllines:
            if 'Cell volume [ang^3]' in i:
                cell_volumeli = i.split()
                cell_volumelis.append(float(cell_volumeli[-1]))
        cell_volumelis=[i*1E-3 for i in cell_volumelis]




        fig,(ax1, ax2) = plt.subplots(1,2,figsize=(17.5*self.inch_cm,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)
        # ax1.scatter(steplis, pressureint, alpha=0.85, marker='o',color='NONE',edgecolor='red', s=5,linewidth=1)
        ax1.plot(steplis, pressureint,alpha=0.75,linewidth=2.8,label='Pressure')
        ax2.plot(steplis,cell_volumelis,alpha=0.75,linewidth=2.4,label='Cell volume')

        # plt.axhline(y=2,c="black")
        # plt.axvline(x=534,c="green")#461-500,520-534,
        # 以下是图片的格式设置
        # 设置横纵坐标的名称以及对应字体格式
        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'8'
        }
        # ax.legend(prop=font1,frameon=False,)
        # ax.set_xlim((320,600))
        ax1.set_xlabel(r'Step',fontsize=8,family='Arial')
        ax1.set_ylabel(r'Pressure (atm)',fontsize=8,family='Arial')
        ax2.set_xlabel(r'Step',fontsize=8,family='Arial')
        ax2.set_ylabel(r'Cell volume (nm $^3$)',fontsize=8,family='Arial')

        # ax1 = ax.twiny()
        # ax1.scatter(times, weight02,marker='o',alpha=0.95,color='NONE', s=5,edgecolors='blue',linewidth=0.8,label='Weight-Time')
        # font1 = {'family' : 'Arial',
        # 'weight' : 'normal','size':'8'
        # }
        ax1.legend(prop=font1,frameon=False,)
        # ax.set_xlim((-10,270))
        # ax1.set_xlabel(r'Time (min)',fontsize=8,family='Arial')
        # ax1.set_ylabel(r'Weight (%)',fontsize=8,family='Arial')
        # cd=max(a5)/10
        # ax.set_xlim((max(wavenumber),min(wavenumber)))
        # ax1.set_xlim((20,40))
        #不显示Y轴的刻度
        ax2.legend(prop=font1,frameon=False,ncol=1)
        ax1.set_xticks([i for i in np.arange(0,len(steplis),len(steplis)//5)])
        ax2.set_xticks([i for i in np.arange(0,len(steplis),len(steplis)//5)])
        # plt.yticks([i for i in np.arange(0,1.01,0.25)])
        # ax1.set_yticks([])
        # ax2.set_yticks([])
        #存储为
        os.chdir(self.dire)
        plt.savefig('pressurechek.tiff', transparent=True,format='tiff')#指定分辨率,边界紧，背景透明bbox_inches='tight',
        plt.show()
        print('\n')
        print('\n')
        print('*******************************************************************************')
        print('\n')
        print('AVERAGE PRESSURE [bar]               =                   %s'%str(np.mean([float(i) for i in pressurelis])))   
        print('FINAL PRESSURE [bar]                 =                   %s'%str(np.mean([float(i) for i in pressurelis][-1])))
        print('\n')
        print('*******************************************************************************')
        print ('congratulate！！！')


    def frame(self,filename='aimd.cif',sccs=0,semethod='QS',runtype='MD'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        import re
        os.chdir(self.dire)
######################################################################################################################
        element=['Ca','H','S','O','Na','Fe','S','O']
        Valence_electron_number=[10,1,6,6,1,2,3,5]
####################################################################################################################
        inputli=['\n','100','2','33','\n','0','q']
        filenas = str('tes.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd='Multiwfn %s -isilent < tes.txt'%str(filename)
        print(cmd)
        res=os.popen(cmd)
        output_str=res.read()
        print(output_str)
        os.remove(filenas)   
######################################################################################################################
        inputline=[]
#SET VARIABLE
        line1_1='!!!   Written by Wcl    !!!'
        inputline.append(line1_1)
        line1_2=' '
        inputline.append(line1_2)
        line1_3='@SET project_name aimd'
        inputline.append(line1_3)
        line1_04='@SET run_type %s'%str(runtype)##########         此处更改任务类型
        inputline.append(line1_04)
        line1_4='@SET FORCE_METHOD %s'%str(semethod)######        此处更改计算方法
        inputline.append(line1_4)
        line1_5='@SET RESTART 0'
        inputline.append(line1_5)
        line1_6='@SET COLVAR 0'
        inputline.append(line1_6)

        if '1' in line1_6:
            print('                 Please be sure to provide "colvar.inp" file !!!              ')
        line1_7='@SET VELOCITY 0'
        inputline.append(line1_7)
        if '1' in line1_7:
            print('                 Please be sure to provide "velocity.inp" file !!!              ')
        line1_8='@SET SCCS %s'%str(sccs)
        inputline.append(line1_8)

        line1__1=' '
        inputline.append(line1__1)        
######################################################################################################################
#GLOBAL
        line2_1='&GLOBAL'
        inputline.append(line2_1)
        line2_2='  PROJECT ${project_name}'
        inputline.append(line2_2)
        line2_3='  PRINT_LEVEL LOW'
        inputline.append(line2_3)
        line2_4='  RUN_TYPE ${run_type}'
        inputline.append(line2_4)
        line2_5='&END GLOBAL'
        inputline.append(line2_5)

######################################################################################################################
#FORCE_EVAL_FRAME
        line3_1='&FORCE_EVAL'
        inputline.append(line3_1)
        line3_2='  STRESS_TENSOR ANALYTICAL'
        inputline.append(line3_2)
        line3_3='  METHOD ${FORCE_METHOD}'
        inputline.append(line3_3)
#       ##########################################################################################################
#       FORCE_EVAL_SUBSYS
        line3_4='  &SUBSYS'
        inputline.append(line3_4)
        line3_5='    &CELL'
        inputline.append(line3_5)
        structurefilename=filename
        structure = Structure.from_file(structurefilename)
        latt=structure.lattice.as_dict()
        latt=latt.get('matrix')
        line3_6='      A    '+'%s     '%str(format(latt[0][0],'10f'))+'%s     '%str(format(latt[0][1],'10f'))+'%s     '%str(format(latt[0][2],'10f'))
        inputline.append(line3_6)
        line3_7='      B    '+'%s     '%str(format(latt[1][0],'10f'))+'%s     '%str(format(latt[1][1],'10f'))+'%s     '%str(format(latt[1][2],'10f'))
        inputline.append(line3_7)
        line3_8='      C    '+'%s     '%str(format(latt[2][0],'10f'))+'%s     '%str(format(latt[2][1],'10f'))+'%s     '%str(format(latt[2][2],'10f'))
        inputline.append(line3_8)
        line3_9='      PERIODIC XYZ'
        inputline.append(line3_9)
        line3_10='    &END CELL'
        inputline.append(line3_10)

        line3_10='    &TOPOLOGY'
        inputline.append(line3_10)
        line3_10='      COORD_FILE_NAME ./aimd.pdb'
        inputline.append(line3_10)
        line3_10='      COORDINATE pdb'
        inputline.append(line3_10)
        line3_10='      CONN_FILE_FORMAT OFF'
        inputline.append(line3_10)        
        line3_10='    &END TOPOLOGY'
        inputline.append(line3_10)
        
        line3_16='    @if ${COLVAR} == 1'
        inputline.append(line3_16)
        line3_17='    @include "colvar.inp"'
        inputline.append(line3_17)    
        line3_18='    @endif'
        inputline.append(line3_18) 


        line3_19='    @if ${VELOCITY} == 1'
        inputline.append(line3_19)
        line3_20='    @include "velocity.inp"'
        inputline.append(line3_20)    
        line3_21='    @endif'
        inputline.append(line3_21) 


        dicts=list(structure.as_dict().items())[4][1]

        # print(dicts[0]['xyz'])
###################################################################################################################
        # Atomic_charge_element={'Ca':2,'H':1,'S':-2,'O':-2,'Na':1,'Fe':2}
        abc=structure.lattice.abc
        angles=structure.lattice.angles
        stdate=structure.as_dataframe()
        extract=stdate.drop_duplicates(subset='Species',keep='first')
        extract=np.array(extract['Species']).tolist()
        extract=[str(i)[:-1] for i in extract]
        atomleble=[]
        elements=[]
        for n, u in enumerate(dicts):
            xx=np.array(u['xyz'])
            dictscal=[dicts[f] for f in np.arange(len(dicts)) if f != n and dicts[f]['label']!=u['label']]
            resus=[]
            for j,g in enumerate(dictscal):
                yy=np.array(g['xyz'])
                resu=np.linalg.norm(xx-yy)
                resus.append(resu)
            inde=resus.index(min(resus))
            if dictscal[inde]['label']=='H':
                atomle=str(u['label']).upper()+str(dictscal[inde]['label']).upper()
            else:
                atomle=str(u['label']).upper()
            element=u['label']
            atomleble.append(atomle)
            elements.append(element)

        # print(stdate['Species'])
        # chargelist=[Atomic_charge_element[str(i)[:-1]] for i in (stdate['Species'])]


        data = {'atom':['ATOM' for i in np.arange(len(stdate['Species']))],
                'index':[i+1 for i in np.arange(len(stdate['Species']))],
                'name':atomleble,
                'x':['x' for i in np.arange(len(stdate['Species']))],
                '1':[1 for i in np.arange(len(stdate['Species']))],
                '11':['1.00' for i in np.arange(len(stdate['Species']))],
                '00':['0.00' for i in np.arange(len(stdate['Species']))],
                'charges':['  ' for i in np.arange(len(stdate['Species']))],
                'element':elements}

        atomfield = pd.DataFrame(data)
        pdbfield=pd.concat([atomfield.iloc[:,0:5],stdate.iloc[:,4:7],atomfield.iloc[:,5:7],atomfield.iloc[:,8],atomfield.iloc[:,7]],axis=1)
        stdate=structure.as_dataframe()
        extract=pdbfield.drop_duplicates(subset='name',keep='first')
        extract_name=np.array(extract['name']).tolist()
        extract_element=np.array(extract['element']).tolist()

#############################################################################################################################################

        if 'QS' in line1_4:
            for i,h in enumerate(extract_name):
                line3_11='    &KIND %s'%str(extract_name[i])
                inputline.append(line3_11)
                line3_12='      ELEMENT %s'%str(extract_element[i])
                inputline.append(line3_12)
                line3_13='      BASIS_SET DZVP-MOLOPT-SR-GTH'
                inputline.append(line3_13)
                line3_14='      POTENTIAL GTH-PBE'
                inputline.append(line3_14)
                line3_15='    &END KIND'
                inputline.append(line3_15)
            line3__1='  &END SUBSYS'
            inputline.append(line3__1)
            line3_15='  @include "dft.inp"'
            inputline.append(line3_15)
            print('                 Please be sure to provide "dft.inp" file !!!              ')

        if 'FIST' in line1_4:
            for i,h in enumerate(extract_name):
                line3_11='    &KIND %s'%str(extract_name[i])
                inputline.append(line3_11)
                line3_12='      ELEMENT %s'%str(extract_element[i])
                inputline.append(line3_12)
                line3_15='    &END KIND'
                inputline.append(line3_15)
            line3__1='  &END SUBSYS'
            inputline.append(line3__1)
            line3_15='  @include "mm.inp"'
            inputline.append(line3_15)
            print('                 Please be sure to provide "mm.inp" file !!!              ')

        if 'QMMM' in line1_4:
            for i,h in enumerate(extract_name):
                line3_11='    &KIND %s'%str(extract_name[i])
                inputline.append(line3_11)
                line3_12='      ELEMENT %s'%str(extract_element[i])
                inputline.append(line3_12)
                line3_13='      BASIS_SET DZVP-MOLOPT-SR-GTH'
                inputline.append(line3_13)
                line3_14='      POTENTIAL GTH-PBE'
                inputline.append(line3_14)
                line3_15='    &END KIND'
                inputline.append(line3_15)
            line3__1='  &END SUBSYS'
            inputline.append(line3__1)
            line3_15='  @include "dft.inp"'
            inputline.append(line3_15)
            print('                 Please be sure to provide "dft.inp" file !!!              ')            
            line3_15='  @include "mm.inp"'
            inputline.append(line3_15)
            print('                 Please be sure to provide "mm.inp" file !!!              ')

        line3_1='&END FORCE_EVAL'
        inputline.append(line3_1)
######################################################################################################################
#MOTION
        if 'MD' in line1_04:
            line4_1='@include "motion.inp"'
            inputline.append(line4_1)
        elif 'OPT' in line1_04:
            line4_1='@include "motion.inp"'
            inputline.append(line4_1)
        elif 'VIBRATIONAL_ANALYSIS' in line1_04:
            line4_1='@include "motion.inp"'
            inputline.append(line4_1)
        else:
            print(format('                 The "motion.inp" field is not required',' <15'))


######################################################################################################################
#RESTART
        line9_0=' '
        inputline.append(line9_0) 
        line9_1='@if ${RESTART} == 1'
        inputline.append(line9_1)
        line9_2='&EXT_RESTART'
        inputline.append(line9_2)  
        line9_3='  RESTART_FILE_NAME ./${project_name}-1.restart'
        inputline.append(line9_3)    
        line9_4='  RESTART_DEFAULT T'
        inputline.append(line9_4)  
        line9_5='&END EXT_RESTART'
        inputline.append(line9_5)  
        line9_6='@endif'
        inputline.append(line9_6)                                
######################################################################################################################
#WRITE
        filename = str('aimd.inp')
        s='\n'.join(inputline)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()


    def pdbgenerate(self,filename='aimd.cif'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        os.chdir(self.dire)
######################################################################################################################
        inputli=['\n','100','2','1','\n','0','q']
        filenas = str('tes.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd='Multiwfn %s -isilent < tes.txt'%str(filename)
        # print(cmd)
        res=os.popen(cmd)
        output_str=res.read()
        # print(output_str)
        os.remove(filenas)  
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')
        filename='aimd.pdb'
######################################################################################################################
        f=open(filename,'r')
        alllines=f.readlines()#.strip().split()
        f.close()
        alls=alllines[2:-1]
        alllines[0] = alllines[0][:9] + 'Written by wcl' + alllines[0][30:]
        for n, u in enumerate(alls):
            uvalue=u.split()
            xx=np.array([float(i) for i in uvalue[-6:-3]])
            dictscal=[alls[f] for f in np.arange(len(alls)) if f != n and alls[f][2]!=uvalue[2]]
            resus=[]
            for j,g in enumerate(dictscal):
                gvalue=g.split()
                yy=np.array([float(i) for i in gvalue[-6:-3]])
                # print('###########%s'%(yy))
                resu=np.linalg.norm(xx-yy)
                resus.append(resu)
            inde=resus.index(min(resus))
            
            if dictscal[inde].split()[2]=='H':
                atomle=(str(uvalue[2]).upper())+str(dictscal[inde].split()[2]).upper()
                atomle=str(format(str(atomle),' <2'))
                # print(atomle)
                alllines[n+2] = alllines[n+2][:13] + atomle + alllines[n+2][15:]           
            else:
                atomle=str(uvalue[2]).upper()
                atomle=str(format(str(atomle),' <2'))
                # print(atomle)
                alllines[n+2] = alllines[n+2][:13] + atomle + alllines[n+2][15:]
            print(alllines[n+2])               
######################################################################################################################
#WRITE
        filename = str('aimd.pdb')
        s=''.join(alllines)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()



    def motiongenerate(self,ensemble='NPT_F',fixatoms=['Ca','S']):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        import re
        os.chdir(self.dire)
        filename='aimd.cif'
######################################################################################################################
        inputli=['\n','100','2','1','\n','0','q']
        filenas = str('tes.txt')
        sg='\n'.join(inputli)
        fp=open(filenas,'w+')
        fp.write(sg)
        fp.close()
######################################################################################################################
#perform Multiwfn
        cmd='Multiwfn %s -isilent < tes.txt'%str(filename)
        print(cmd)
        res=os.popen(cmd)
        output_str=res.read()
        print(output_str)
        os.remove(filenas)  
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')
        filename='aimd.pdb'
######################################################################################################################
        f=open(filename,'r')
        alllines=f.readlines()#.strip().split()
        f.close()
        alls=alllines[2:-1]
        alllines[0] = alllines[0][:9] + 'Written by wcl' + alllines[0][30:]
        atomindexfix=[]
        for n, u in enumerate(alls):
            uvalue=u.split()
            xx=np.array([float(i) for i in uvalue[-6:-3]])
            dictscal=[alls[f] for f in np.arange(len(alls)) if f != n and alls[f][2]!=uvalue[2]]
            resus=[]
            for j,g in enumerate(dictscal):
                gvalue=g.split()
                yy=np.array([float(i) for i in gvalue[-6:-3]])

                resu=np.linalg.norm(xx-yy)
                resus.append(resu)
            inde=resus.index(min(resus))
            if type(fixatoms) == list:
                jhg=1
                # print(u)
                if uvalue[2] in '%s'%str(fixatoms):# or u['label']=='%s'%str(fixatoms[1]) or u['label']=='%s'%str(fixatoms[2]):
                    # print(uvalue[2])
                    atomindexfix.append(str(n+1))                
            elif type(fixatoms) == tuple:
                jhg=1
                if xx[2]<fixatoms[0] or xx[2]>fixatoms[1]:
                    atomindexfix.append(str(n+1))
            else:
                jhg=0
        if jhg==0:
            print("\033[1;30;47m\tPlease define your own fixed atom\033[0m")
        ss=' '.join(atomindexfix)
        print('\033[1;30;47m\t\n\n# The fixed list of atoms is #\n %s \n\033[0m'%str(ss))
######################################################################################################################
        inputline=[]
        line1_1='!!!   Written by Wcl    !!!'
        inputline.append(line1_1)
        line1_2=' '
        inputline.append(line1_2)
        line1_4='@SET Ensemble %s'%str(ensemble)##########         此处更改任务类型
        inputline.append(line1_4)
        line1_2=' '
        inputline.append(line1_2)    
######################################################################################################################
#       MOTION FOR MD 


        line2_1='&MOTION'
        inputline.append(line2_1)
        line2_2='  &CONSTRAINT'
        inputline.append(line2_2)
        line2_3='     &FIXED_ATOMS'
        inputline.append(line2_3)
        line2_4='     LIST %s'%str(ss)
        inputline.append(line2_4)
        line2_5='     &END FIXED_ATOMS'
        inputline.append(line2_5)
        line2_2='  &END CONSTRAINT'
        inputline.append(line2_2)
        line3_1='  &PRINT'
        inputline.append(line3_1)
        line3_2='    &FORCES'
        inputline.append(line3_2)
        line3_3='      ADD_LAST NUMERIC'
        inputline.append(line3_3)
        line3_4='    &END FORCES'
        inputline.append(line3_4)
        line3_5='    &TRAJECTORY'
        inputline.append(line3_5)
        line3_6='          LOG_PRINT_KEY T'
        inputline.append(line3_6)
        line3_7='          FORMAT PDB'
        inputline.append(line3_7)
        line3_4='          &EACH'
        inputline.append(line3_4)
        line3_5='            MD 1'
        inputline.append(line3_5)
        line3_6='          &END EACH'
        inputline.append(line3_6)
        line3_7='          ADD_LAST NUMERIC'
        inputline.append(line3_7)

        line3_4='    &END TRAJECTORY'
        inputline.append(line3_4)
        line3_5='    &CELL'
        inputline.append(line3_5)
        line3_6='    ADD_LAST NUMERIC'
        inputline.append(line3_6)
        line3_7='    &END CELL'
        inputline.append(line3_7)

        line3_2='    &VELOCITIES'
        inputline.append(line3_2)
        line3_3='      &EACH'
        inputline.append(line3_3)
        line3_4='        MD     1 '
        inputline.append(line3_4)
        line3_5='      &END EACH'
        inputline.append(line3_5)
        line3_6='    &END VELOCITIES'
        inputline.append(line3_6)
        line3_7='    &RESTART'
        inputline.append(line3_7)
        line3_4='      BACKUP_COPIES 0'
        inputline.append(line3_4)
        line3_5='      &EACH'
        inputline.append(line3_5)
        line3_6='        MD 10'
        inputline.append(line3_6)
        line3_7='      &END EACH'
        inputline.append(line3_7)
        line3_4='    &END RESTART'
        inputline.append(line3_4)
        line3_7='    &RESTART_HISTORY'
        inputline.append(line3_7)
        line3_5='      &EACH'
        inputline.append(line3_5)
        line3_6='        MD 0'
        inputline.append(line3_6)
        line3_6='        CELL_OPT 0'
        inputline.append(line3_6)
        line3_7='      &END EACH'
        inputline.append(line3_7)
        line3_4='    &END RESTART_HISTORY'
        inputline.append(line3_4)

        line3_5='  &END PRINT'
        inputline.append(line3_5)

        line3_16='@if ${COLVAR} == 1'
        inputline.append(line3_16)
        line3_1='  &FREE_ENERGY '
        inputline.append(line3_1)
        line3_2='    &METADYN'
        inputline.append(line3_2)
        line3_3='      DO_HILLS .FALSE.'
        inputline.append(line3_3)
        line3_4='      &METAVAR'
        inputline.append(line3_4)
        line3_5='        SCALE 0.2'
        inputline.append(line3_5)
        line3_6='        COLVAR 1'
        inputline.append(line3_6)
        line3_7='      &END METAVAR'
        inputline.append(line3_7)

        line3_4='      &METAVAR'
        inputline.append(line3_4)
        line3_5='        SCALE 0.2'
        inputline.append(line3_5)
        line3_6='        COLVAR 2'
        inputline.append(line3_6)
        line3_7='      &END METAVAR'
        inputline.append(line3_7)

        line3_4='      &PRINT'
        inputline.append(line3_4)
        line3_5='        &COLVAR'
        inputline.append(line3_5)
        line3_6='          COMMON_ITERATION_LEVELS 3'
        inputline.append(line3_6)
        line3_7='          &EACH'
        inputline.append(line3_7)
        line3_4='            MD 1'
        inputline.append(line3_4)
        line3_5='          &END'
        inputline.append(line3_5)
        line3_6='        &END COLVAR'
        inputline.append(line3_6)
        line3_7='      &END PRINT'
        inputline.append(line3_7)
        line3_2='    &END METADYN'
        inputline.append(line3_2)
        line3_1='  &END FREE_ENERGY '
        inputline.append(line3_1)
        line3_18='@endif'
        inputline.append(line3_18)

        line2_1='@if ${run_type} == MD'
        inputline.append(line2_1)
        line2_2='  &MD'
        inputline.append(line2_2)
        line2_3='    ENSEMBLE ${Ensemble}'
        inputline.append(line2_3)
        line2_4='    STEPS 1000'
        inputline.append(line2_4)
        line2_5='    TIMESTEP 0.5'
        inputline.append(line2_5)

        line3_1='    TEMPERATURE 298.15'
        inputline.append(line3_1)
        line3_2='    &THERMOSTAT'
        inputline.append(line3_2)
        line3_3='      TYPE CSVR'
        inputline.append(line3_3)
        line3_4='      &CSVR'
        inputline.append(line3_4)
        line3_5='        TIMECON 300'
        inputline.append(line3_5)
        line3_6='      &END CSVR'
        inputline.append(line3_6)
        line3_7='    &END THERMOSTAT'
        inputline.append(line3_7)

        if 'NPT' in line1_4:
            line3_2='    &BAROSTAT'
            inputline.append(line3_2)
            line3_3='      PRESSURE 1.01325'
            inputline.append(line3_3)
            line3_4='      VIRIAL XYZ'
            inputline.append(line3_4)
            line3_5='    &END BAROSTAT'
            inputline.append(line3_5)
            print('                 Now the pressure under the NPT ensemble is 1.01325               ')
        line3_8='  &END MD'
        inputline.append(line3_8)
        line2_1='@endif'
        inputline.append(line2_1)     

######################################################################################################################
#       MOTION FOR OPT 
        line1_2=' '
        inputline.append(line1_2) 
        line2_1='@if ${run_type} == CELL_OPT'
        inputline.append(line2_1)
        line2_2='  &CELL_OPT'
        inputline.append(line2_2)
        line2_3='    EXTERNAL_PRESSURE 1.01325 #External pressure for cell optimization (bar)'
        inputline.append(line2_3)
        line2_4='    CONSTRAINT NONE #Can be e.g. Z, XY to fix corresponding cell length'
        inputline.append(line2_4)
        line2_5='    KEEP_ANGLES F #If T, then cell angles will be kepted'
        inputline.append(line2_5)

        line3_1='    KEEP_SYMMETRY F #If T, crystal symmetry will be kepted, and symmetry should be specified in &CELL'
        inputline.append(line3_1)
        line3_2='    TYPE DIRECT_CELL_OPT #Geometry and cell are optimized at the same time. Can also be GEO_OPT, MD'
        inputline.append(line3_2)
        line3_3='    MAX_DR 3E-3 #Maximum geometry change'
        inputline.append(line3_3)
        line3_4='    RMS_DR 1.5E-3 #RMS geometry change'
        inputline.append(line3_4)
        line3_5='    MAX_FORCE 4.5E-4 #Maximum force'
        inputline.append(line3_5)
        line3_6='    RMS_FORCE 3E-4 #RMS force'
        inputline.append(line3_6)
        line3_7='    PRESSURE_TOLERANCE 100 #Pressure tolerance (w.r.t EXTERNAL_PRESSURE)'
        inputline.append(line3_7)


        line3_2='    OPTIMIZER CG #Can also be BFGS or LBFGS'
        inputline.append(line3_2)
        line3_3='    &CG'
        inputline.append(line3_3)
        line3_4='      &LINE_SEARCH'
        inputline.append(line3_4)
        line3_5='        TYPE 2PNT #Two-point extrapolation, cheap while acceptable. Can also be FIT, GOLD'
        inputline.append(line3_5)
        line3_8='      &END LINE_SEARCH'
        inputline.append(line3_8)
        line3_2='    &END CG'
        inputline.append(line3_2)
        line3_3='  &END CELL_OPT'
        inputline.append(line3_3)
        line2_1='@endif'
        inputline.append(line2_1)
        line3_6='&END MOTION'
        inputline.append(line3_6)

######################################################################################################################
#       MOTION FOR VIBRATIONAL_ANALYSIS 
        line1_2=' '
        inputline.append(line1_2) 

        line2_1='@if ${run_type} == VIBRATIONAL_ANALYSIS'
        inputline.append(line2_1)
        line3_6='&VIBRATIONAL_ANALYSIS'
        inputline.append(line3_6)
        line2_2='  DX 0.01 #Step size of finite difference. This is default (Bohr)'
        inputline.append(line2_2)
        line2_3='  NPROC_REP 64 #Number of processors to be used per replica. This is default'
        inputline.append(line2_3)
        line2_4='  TC_PRESSURE 101325 #1 atm. Pressure for calculate thermodynamic data (Pa)'
        inputline.append(line2_4)
        line2_5='  TC_TEMPERATURE 298.15 #Temperature for calculate thermodynamic data (K)'
        inputline.append(line2_5)

        line3_1='  THERMOCHEMISTRY #Print thermochemistry information (only correct for gas molecule!)'
        inputline.append(line3_1)
        line3_2='  INTENSITIES T #Calculate IR intensities'
        inputline.append(line3_2)
        '''
        line3_3='  &MODE_SELECTIVE'
        inputline.append(line3_3)
        line3_4='    ATOMS %s %s #Specifies the list of atoms which should be displaced for the Initial guess'%(fixatoms[0],fixatoms[1])
        inputline.append(line3_4)    
        line3_5='      INVOLVED_ATOMS %s %s'%(fixatoms[0],fixatoms[1])
        inputline.append(line3_5)
        line3_5='    &END INVOLVED_ATOMS'
        inputline.append(line3_5)
        line3_6='  &END MODE_SELECTIVE'
        inputline.append(line3_6)
        '''
        line3_3='  &PRINT'
        inputline.append(line3_3)
        line3_4='    &MOLDEN_VIB #Output .mol (Molden file) for visualization vibrational modes'
        inputline.append(line3_4)
        line3_5='    &END MOLDEN_VIB'
        inputline.append(line3_5)
        line3_6='  &END PRINT'
        inputline.append(line3_6)
        line3_7='&END VIBRATIONAL_ANALYSIS'
        inputline.append(line3_7)
        line2_1='@endif'
        inputline.append(line2_1)
                                
######################################################################################################################
#WRITE
        wfilename = str('motion.inp')
        s='\n'.join(inputline)
        fp=open(wfilename,'w+')
        fp.write(s)
        fp.close()


    def dftgenerate(self,method='OT',ldoslist=[1,2,3]):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        import re
        os.chdir(self.dire)
######################################################################################################################
        element=['Ca','H','S','O']
        Valence_electron_number=[10,1,6,6]


######################################################################################################################
        inputline=[]
        line1_1='!!!   Written by Wcl    !!!'
        inputline.append(line1_1)
        line1_2=' '
        inputline.append(line1_2)  
######################################################################################################################
#       MOTION FOR MD 


        line2_1='&DFT'
        inputline.append(line2_1)

        line2_1='  BASIS_SET_FILE_NAME  BASIS_MOLOPT'
        inputline.append(line2_1)
        line2_2='  POTENTIAL_FILE_NAME  POTENTIAL'
        inputline.append(line2_2)
        line2_3='#  WFN_RESTART_FILE_NAME aimd-RESTART.wfn'
        inputline.append(line2_3)
        line2_4='  CHARGE    0 #Net charge'
        inputline.append(line2_4)
        line2_5='  MULTIPLICITY    1 '
        inputline.append(line2_5)

        line3_1='  &QS'
        inputline.append(line3_1)
        line3_2='  EPS_DEFAULT 1E-10 #This is default. Set all EPS_xxx to values such that the energy will be correct up to this value'
        inputline.append(line3_2)
        line3_3='  &END QS'
        inputline.append(line3_3)
        line3_4='  &POISSON'
        inputline.append(line3_4)
        line3_5='    PERIODIC XYZ #Direction(s) of PBC for calculating electrostatics'
        inputline.append(line3_5)
        line3_6='    PSOLVER PERIODIC #The way to solve Poisson equation'
        inputline.append(line3_6)
        line3_7='  &END POISSON'
        inputline.append(line3_7)

        line3_2='  &XC'
        inputline.append(line3_2)
        line3_3='    &XC_FUNCTIONAL PBE'
        inputline.append(line3_3)
        line3_4='    &END XC_FUNCTIONAL'
        inputline.append(line3_4)
        line3_5='    &VDW_POTENTIAL'
        inputline.append(line3_5)
        line3_8='      POTENTIAL_TYPE PAIR_POTENTIAL'
        inputline.append(line3_8)
        line3_16='     &PAIR_POTENTIAL'
        inputline.append(line3_16)
        line3_1='        PARAMETER_FILE_NAME dftd3.dat'
        inputline.append(line3_1)
        line3_2='        TYPE DFTD3'
        inputline.append(line3_2)
        line3_3='        REFERENCE_FUNCTIONAL PBE'
        inputline.append(line3_3)
        line3_4='      &END PAIR_POTENTIAL'
        inputline.append(line3_4)
        line3_5='    &END VDW_POTENTIAL'
        inputline.append(line3_5)
        line3_6='  &END XC'
        inputline.append(line3_6)
        line3_7='  &MGRID'
        inputline.append(line3_7)
        line3_4='    CUTOFF 400'
        inputline.append(line3_4)
        line3_5='    REL_CUTOFF 55'
        inputline.append(line3_5)
        line3_6='  &END MGRID'
        inputline.append(line3_6)

        if "OT" in method:
            line3_7='  &SCF'
            inputline.append(line3_7)
            line3_4='    MAX_SCF 20 #Maximum number of steps of inner SCF'
            inputline.append(line3_4)
            line3_5='    EPS_SCF 3.0E-06 #Convergence threshold of density matrix of inner SCF'
            inputline.append(line3_5)
            line3_6='#   SCF_GUESS RESTART #Use wavefunction from WFN_RESTART_FILE_NAME file as initial guess'
            inputline.append(line3_6)
            line3_7='    &OT'
            inputline.append(line3_7)
            line3_4='      PRECONDITIONER FULL_KINETIC'
            inputline.append(line3_4)
            line3_5='      MINIMIZER DIIS '
            inputline.append(line3_5)
            line3_6='      LINESEARCH 2PNT '
            inputline.append(line3_6)
            line3_7='    &END OT'
            inputline.append(line3_7)
            line3_2='     &OUTER_SCF'
            inputline.append(line3_2)
            line3_1='       MAX_SCF 20 '
            inputline.append(line3_1)
            line3_18='      EPS_SCF 3.0E-06 '
            inputline.append(line3_18)
            line3_1='     &END OUTER_SCF'
            inputline.append(line3_1)
            line2_1='  &END SCF'
            inputline.append(line2_1)
        elif "DI" in method:
            line3_7='  #&KPOINTS'
            inputline.append(line3_7)
            line3_4='    #SCHEME MONKHORST-PACK  1  1  1'
            inputline.append(line3_4)
            line3_4='    #SYMMETRY F #If using symmetry to reduce the number of k-points'
            inputline.append(line3_4)
            line3_7='  #&END KPOINTS'
            inputline.append(line3_7)  
            line3_7='  #UKS'
            inputline.append(line3_7)                       
            line3_7='  &SCF'
            inputline.append(line3_7)
            line3_4='    MAX_SCF 128 #Maximum number of steps of inner SCF'
            inputline.append(line3_4)
            line3_5='    EPS_SCF 3.0E-06 #Convergence threshold of density matrix of inner SCF'
            inputline.append(line3_5)
            line3_6='#   SCF_GUESS RESTART #Use wavefunction from WFN_RESTART_FILE_NAME file as initial guess'
            inputline.append(line3_6)
            line3_7='    &DIAGONALIZATION'
            inputline.append(line3_7)
            line3_4='      ALGORITHM STANDARD'
            inputline.append(line3_4)
            line3_7='    &END DIAGONALIZATION'
            inputline.append(line3_7)
            line3_2='    &MIXING'
            inputline.append(line3_2)
            line3_1='      METHOD BROYDEN_MIXING #PULAY_MIXING is also a good alternative'
            inputline.append(line3_1)
            line3_18='      ALPHA 0.1 #Default. Mixing forty percent of new density matrix with the old one'
            inputline.append(line3_18)
            line3_18='      NBROYDEN 8 #Default is 4. Number of previous steps stored for the actual mixing scheme'
            inputline.append(line3_18)
            line3_1='    &END MIXING'
            inputline.append(line3_1)
            line3_2='    &SMEAR'
            inputline.append(line3_2)
            line3_1='      METHOD FERMI_DIRAC'
            inputline.append(line3_1)
            line3_18='      ELECTRONIC_TEMPERATURE 300 #Electronic temperature of Fermi-Dirac smearing in K'
            inputline.append(line3_18)
            line3_1='    &END SMEAR'
            inputline.append(line3_1)
            line3_1='    ADDED_MOS 30 #Number of virtual MOs to be solved'
            inputline.append(line3_1)
            line2_1='  &END SCF'
            inputline.append(line2_1)            

        line3_16=' @if ${run_type} == MD'
        inputline.append(line3_16)
        line3_2='  &PRINT'
        inputline.append(line3_2)
        line3_3='     &MULLIKEN'
        inputline.append(line3_3)
        line3_4='         &EACH'
        inputline.append(line3_4)
        line3_4='             MD  1'
        inputline.append(line3_4)        
        line3_4='         &END EACH'
        inputline.append(line3_4)
        line3_4='     &END MULLIKEN'
        inputline.append(line3_4)
        line3_3='     &HIRSHFELD'
        inputline.append(line3_3)
        line3_4='         SHAPE_FUNCTION DENSITY'
        inputline.append(line3_4)        
        line3_4='         &EACH'
        inputline.append(line3_4)
        line3_4='             MD  1'
        inputline.append(line3_4)        
        line3_4='         &END EACH'
        inputline.append(line3_4)
        line3_4='     &END HIRSHFELD'
        inputline.append(line3_4)

        line3_7='  &END PRINT'
        inputline.append(line3_7)
        line3_18=' @endif'
        inputline.append(line3_18)         



        line3_16=' @if ${SCCS} == 1'
        inputline.append(line3_16)
        line3_2='  &SCCS'
        inputline.append(line3_2)

        line3_3='    ALPHA [N*m^-1] 0.0'
        inputline.append(line3_3)
        line3_4='    BETA [kbar] 0.0'
        inputline.append(line3_4)
        line3_5='    GAMMA [mN/m] 0.0'
        inputline.append(line3_5)
        line3_6='    DIELECTRIC_CONSTANT 78.36'
        inputline.append(line3_6)

        line3_3='    EPS_SCCS 1E-4'
        inputline.append(line3_3)
        line3_4='    EPS_SCF 0.5'
        inputline.append(line3_4)
        line3_3='    MAX_ITER 100'
        inputline.append(line3_3)
        line3_7='    DERIVATIVE_METHOD CD3 #Default. Method for calculation of numerical derivatives. Can also be CD3, CD5, CD7'
        inputline.append(line3_7)
        line3_6='    &ANDREUSSI'
        inputline.append(line3_6)

        line3_3='      RHO_MAX 0.00035  #Small values can help convergence to some extent'
        inputline.append(line3_3)
        line3_4='      RHO_MIN 0.00001  #Small values can help convergence to some extent'
        inputline.append(line3_4)
        line3_3='    &END ANDREUSSI'
        inputline.append(line3_3)
        line3_7='  &END SCCS'
        inputline.append(line3_7)
        line3_18=' @endif'
        inputline.append(line3_18) 

        line3_16=' @if ${run_type} == VIBRATIONAL_ANALYSIS'
        inputline.append(line3_16)
        line3_2='  &PRINT'
        inputline.append(line3_2)

        line3_3='     &MOMENTS'
        inputline.append(line3_3)
        line3_4='     &END MOMENTS'
        inputline.append(line3_4)
        line3_7='  &END PRINT'
        inputline.append(line3_7)
        line3_18=' @endif'
        inputline.append(line3_18) 

        line3_16=' @if ${run_type} == ENERGY'
        inputline.append(line3_16)
        line3_2='  &PRINT'
        inputline.append(line3_2)

        line3_3='     &MO_MOLDEN'
        inputline.append(line3_3)
        line3_4='       NDIGITS 8'
        inputline.append(line3_4)
        line3_5='       GTO_KIND SPHERICAL'
        inputline.append(line3_5)
        line3_6='     &END MO_MOLDEN'
        inputline.append(line3_6)

        line3_3='     &PDOS'
        inputline.append(line3_3)
        line3_4='       FILENAME ./aimd_pdos'
        inputline.append(line3_4)

        for i in ldoslist:
            line3_5='       &LDOS'
            inputline.append(line3_5)
            line3_6='         LIST %s'%str(i)
            inputline.append(line3_6)
            line3_5='       &END LDOS'
            inputline.append(line3_5)

        line3_3='     &END PDOS'
        inputline.append(line3_3)
        line3_3='     &V_HARTREE_CUBE'
        inputline.append(line3_3)
        line3_3='     &END V_HARTREE_CUBE'
        inputline.append(line3_3)
        line3_3='     &E_DENSITY_CUBE'
        inputline.append(line3_3)
        line3_3='     &END E_DENSITY_CUBE'
        inputline.append(line3_3)

        line3_7='  &END PRINT'
        inputline.append(line3_7)
        line3_18=' @endif'
        inputline.append(line3_18) 

       
        line2_1='&END DFT'
        inputline.append(line2_1)                               
######################################################################################################################
#WRITE
        filename = str('dft.inp')
        s='\n'.join(inputline)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()


    def mmgenerate(self,):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        from itertools import product, combinations, permutations
        import re
        os.chdir(self.dire)
######################################################################################################################
#        bond  matrix
        atomlabels=['CA', 'H', 'S', 'O', 'OH']

        bendpremas=product(atomlabels,atomlabels,atomlabels)
        bendma=[]
        for i in bendpremas:
            if i[0]=='H' and i[1]=='OH' and i[2]=='H':
                newval=[i[0],i[1],i[2],'3.29136','106.597']
            elif i[0]=='O' and i[1]=='S' and i[2]=='O':
                newval=[i[0],i[1],i[2],'15.0','109.47']
            else:
                newval=[i[0],i[1],i[2],'0.0','0.0'] 
            bendma.append(newval)

        # print(list(bendma))






        # bendma=[['H','OH','H','3.29136','113.24'],
        #        ['O','S','O','15.0','109.47'],] #['atom1','atom2','atom3','k','theta0']

        bondpremas=product(atomlabels,atomlabels)
        bondma=[]
        for i in bondpremas:
            if i[0]=='H' and i[1]=='OH':
                newval=[i[0],i[1],'900','0.962']
            elif i[0]=='S' and i[1]=='O':
                newval=[i[0],i[1],'5.0 1.2','1.505']
            else:
                newval=[i[0],i[1],'0.0','0.0'] 
            bondma.append(newval)
        # print(list(bondma))

        chargema=[['CA','2.0'],
                 ['H','0.41'],
                 ['S','1.36'],
                 ['OH','-0.82'],
                 ['O','-0.84']]

        ljma=[['OH','OH','0.0067400','3.165492'],
               ['OH','CA','0.0009500','3.350000']]

        # buckpremas=permutations(atomlabels,2)
        # buckma=[]
        # for i in buckpremas:
        #     if i[0]=='O' and i[1]=='O':
        #         newval=[i[0],i[1],'12534.455133','5.0']
        #     elif i[0]=='OH' and i[1]=='O':
        #         newval=[i[0],i[1],'12534.455133','4.065']
        #     elif i[0]=='O' and i[1]=='CA':
        #         newval=[i[0],i[1],'1815.6986','3.5289']
        #     else:
        #         newval=[i[0],i[1],'0.0','0.0'] 
        #     buckma.append(newval)
        # print(buckma)
        buckma=[['O','O','12534.455133','5.0'],
               ['OH','O','12534.455133','4.065'],
               ['O','CA','1815.6986','3.5289'],
               ['CA','CA','0.0','0.0'],
               ['H','CA','0.0','0.0'],
               ['S','CA','0.0','0.0'],
               ['H','H','0.0','0.0'],
               ['S','S','0.0','0.0'],
               ['H','S','0.0','0.0'],
               ['OH','S','0.0','0.0'],
               ['O','H','0.0','0.0'],
               ['OH','H','0.0','0.0'],
               ['O','S','0.0','0.0']]


#############################################################################################################################################

        inputline=[]
        line1_1='!!!   Written by Wcl    !!!'
        inputline.append(line1_1)
        line1_2=' '
        inputline.append(line1_2)  
######################################################################################################################
#       MM FOR FORCEFIELD 


        line2_1='&MM'
        inputline.append(line2_1)
        line2_1='  &FORCEFIELD'
        inputline.append(line2_1)

        for i in bendma:
            line2_2='    &BEND'
            inputline.append(line2_2)
            line2_3='      KIND HARMONIC'
            inputline.append(line2_3)
            line2_4='      ATOMS %s'%str(str(i[0])+' '+str(i[1])+' '+ str(i[2]))
            inputline.append(line2_4)
            line2_5='      K [rad^-2*eV] %s'%str(i[3])
            inputline.append(line2_5)
            line3_1='      THETA0 [deg] %s'%str(i[4])
            inputline.append(line3_1)
            line3_2='    &END BEND'
            inputline.append(line3_2)

        for i in bondma:
            if ' ' in  i[2]:
                line2_2='    &BOND'
                inputline.append(line2_2)
                line2_3='      KIND MORSE'
                inputline.append(line2_3)
                line2_4='      ATOMS %s'%str(str(i[0])+' '+str(i[1]))
                inputline.append(line2_4)
                line2_5='      K [rad^-2kcalmol] %s'%str(i[2])
                inputline.append(line2_5)
                line3_1='      R0 [angstrom] %s'%str(i[3])
                inputline.append(line3_1)
                line3_2='    &END BOND'
                inputline.append(line3_2)
            else:

                line2_2='    &BOND'
                inputline.append(line2_2)
                line2_3='      KIND HARMONIC'
                inputline.append(line2_3)
                line2_4='      ATOMS %s'%str(str(i[0])+' '+str(i[1]))
                inputline.append(line2_4)
                line2_5='      K [angstrom^-2*eV] %s'%str(i[2])
                inputline.append(line2_5)
                line3_1='      R0 [angstrom] %s'%str(i[3])
                inputline.append(line3_1)
                line3_2='    &END BOND'
                inputline.append(line3_2)


        for i in chargema:
            line2_2='    &CHARGE'
            inputline.append(line2_2)
            line2_4='      ATOM %s'%str(str(i[0]))
            inputline.append(line2_4)
            line2_5='      CHARGE %s'%str(i[1])
            inputline.append(line2_5)
            line3_2='    &END CHARGE'
            inputline.append(line3_2)
        line2_1='    &NONBONDED'
        inputline.append(line2_1)

        for i in ljma:
            line2_2='      &LENNARD-JONES'
            inputline.append(line2_2)
            line2_4='        ATOMS %s'%str(str(i[0])+' '+str(i[1]))
            inputline.append(line2_4)
            line2_3='        EPSILON [eV] %s'%str(i[2])
            inputline.append(line2_3)
            line2_5='        SIGMA [angstrom] %s'%str(i[3])
            inputline.append(line2_5)
            line3_1='        RCUT 9.0'
            inputline.append(line3_1)
            line3_2='      &END LENNARD-JONES'
            inputline.append(line3_2)

        for i in buckma:
            line2_2='      &BUCK4RANGES'
            inputline.append(line2_2)
            line2_4='        ATOMS %s'%str(str(i[0])+' '+str(i[1]))
            inputline.append(line2_4)
            line2_3='        A [eV] %s'%str(i[2])
            inputline.append(line2_3)
            line2_5='        B [angstrom] %s'%str(i[3])
            inputline.append(line2_5)
            line2_5='        C 0.0'
            inputline.append(line2_5)
            line3_1='        R1 [angstrom] 19.0'
            inputline.append(line3_1)
            line3_1='        R2 [angstrom] 29.0'
            inputline.append(line3_1)
            line3_1='        R3 [angstrom] 99.0'
            inputline.append(line3_1)
            line3_2='      &END BUCK4RANGES'
            inputline.append(line3_2)

        line2_1='    &END NONBONDED'
        inputline.append(line2_1)
        line2_1='    &SPLINE'
        inputline.append(line2_1)
        line2_1='      EMAX_SPLINE 50E30'
        inputline.append(line2_1)
        line2_1='    &END SPLINE'
        inputline.append(line2_1)
        line2_1='  &END FORCEFIELD'
        inputline.append(line2_1)

######################################################################################################################
#       MM FOR OTHERS 
        line2_1='  &POISSON'
        inputline.append(line2_1)
        line2_1='    PERIODIC XYZ'
        inputline.append(line2_1)
        line2_1='    POISSON_SOLVER PERIODIC'
        inputline.append(line2_1)
        line2_1='    &EWALD'
        inputline.append(line2_1)
        line2_1='      EWALD_TYPE spme'
        inputline.append(line2_1)
        line2_1='      ALPHA 0.44'
        inputline.append(line2_1)
        line2_1='      GMAX 64'
        inputline.append(line2_1)
        # line2_1='      O_SPLINE 6'
        # inputline.append(line2_1)
        line2_1='    &END EWALD'
        inputline.append(line2_1)
        line2_1='  &END POISSON'
        inputline.append(line2_1)
        line2_1='&END MM'
        inputline.append(line2_1)       

        print ('\n\n#####################################################')
        print ('Please find your force field matrices in the source code')
        print ('#####################################################\n\n')

######################################################################################################################
#WRITE
        filename = str('mm.inp')
        s='\n'.join(inputline)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()
######################################################################################################################
#       MM END


    def colvargenerate(self,atomlist=[[1,3],[2,4]]):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        import re
        os.chdir(self.dire)
######################################################################################################################

        inputline=[]
        line1_1='!!!   Written by Wcl    !!!'
        inputline.append(line1_1)
        line1_2=' '
        inputline.append(line1_2)  
######################################################################################################################
#       MOTION FOR MD 

        for i in atomlist:
            fields=str(i[0])+' '+str(i[1])
            line2_1='&COLVAR '
            inputline.append(line2_1)
            line2_1='  &DISTANCE'
            inputline.append(line2_1)
            line2_2='    ATOMS %s'%str(fields)
            inputline.append(line2_2)
            line2_3='  &END DISTANCE'
            inputline.append(line2_3)
            line2_4='&END COLVAR' 
            inputline.append(line2_4)                               
######################################################################################################################
#WRITE
        filename = str('colvar.inp')
        s='\n'.join(inputline)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()


    def shgenerate(self,Ncores=64):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        import re
        os.chdir(self.dire)
######################################################################################################################
        Number_of_cores=Ncores
        Number_of_nodes=Ncores//64


######################################################################################################################
        inputline=[]
        line1_1='#!/bin/bash'
        inputline.append(line1_1)
        line1_2='#SBATCH -p amd_256 '
        inputline.append(line1_2)  
        line2_1='#SBATCH -N %s'%str(Number_of_nodes)
        inputline.append(line2_1)
        line2_1='#SBATCH -n %s'%str(Number_of_cores)
        inputline.append(line2_1)
        line2_2='#SBATCH -J wcl'
        inputline.append(line2_2)
        line2_3='source /public3/soft/modules/module.sh'
        inputline.append(line2_3)
        line2_4='module load cp2k/8.1-public3'
        inputline.append(line2_4)
        line2_5='srun -n %s cp2k.popt -i aimd.inp'%str(Number_of_cores)
        inputline.append(line2_5)

        print('                 The number of cores is %s'%str(Number_of_cores))

######################################################################################################################
#WRITE
        filename = str('sub.sh')

        s='\n'.join(inputline)
        fp=open(filename,'wb')
        fp.write(s.encode('utf-8'))
        fp.close()

######################################################################################################################
#批处理脚本 
        inputline=[]
        line1_1='#!/bin/bash'
        inputline.append(line1_1)
        line1_2='##This script can submit all tasks under the same level directory '
        inputline.append(line1_2)
        line2_1='cd ..'
        inputline.append(line2_1)  
        line2_1='for i in `ls -F |grep /`'
        inputline.append(line2_1)
        line2_1='do'
        inputline.append(line2_1)
        line2_2='cd $i'
        inputline.append(line2_2)
        line2_3='sbatch sub.sh'
        inputline.append(line2_3)
        line2_4='cd ..'
        inputline.append(line2_4)
        line2_5='done'
        inputline.append(line2_5)

######################################################################################################################
#WRITE
        filename = str('batch.sh')

        s='\n'.join(inputline)
        fp=open(filename,'wb')
        fp.write(s.encode('utf-8'))
        fp.close()

    def phono_pre(self,filename='aimd.cif'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        os.chdir(self.dire)
######################################################################################################################
# read phonopy
        f=open('aimd-1.restart','r')
        alllines=f.readlines()#.strip().split()
        f.close()
        
        indexs=[]
        for i,j  in enumerate(alllines):
            if '&CONSTRAINT' in j:
                indexs.append(i)
            elif '&END CONSTRAINT' in j:
                indexs.append(i)
        if indexs:
            del alllines[indexs[0]:indexs[1]+1]
        indexs=[]
        for i,j  in enumerate(alllines):  
            if '&TOPOLOGY' in j:
                indexs.append(i)
            elif '&END TOPOLOGY' in j:
                indexs.append(i)
        if indexs:
            del alllines[indexs[0]:indexs[1]+1]
        indexs=[]
        for i,j  in enumerate(alllines):  
            if '&KIND OH' in j:
                indexs.append(i)
        if indexs:
            del alllines[indexs[0]:indexs[1]+1]

        inputline=[]
        for i,j  in enumerate(alllines):
            if 'CA ' in j or '&KIND CA' in j:
                j=j.replace("CA", "Ca")
                inputline.append(j)
            else:
                inputline.append(j)



######################################################################################################################
#WRITE
        filename = str('aimd.inp')
        s=''.join(inputline)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()

######################################################################################################################
# perform phonopy
        cmd='phonopy --cp2k -d --dim="1 1 1" -c aimd.inp -v'
        print(cmd)
        res=os.popen(cmd)
        output_str=res.read()
        print(output_str)  
        print('###################################################################\n\n')
        print('                    End of external call                                 ')
        print('\n\n###################################################################')

  
######################################################################################################################
# create dir
        filename='supercell-'
        newfilelis=[]
        filelis=[]
        dirlis=[]
        for root,dirs,files in os.walk(os.getcwd()):
            filelis.append(files)
            dirlis.append(dirs)
        filelis=filelis[0]
        dirlis=dirlis[0]


        for i in filelis:
            if filename in i:
                if filename not in str(dirlis):
                    cmd_01='md %s'%(str(i[:-4]))
                    res=os.popen(cmd_01)
                    output_str=res.read()
                    print(output_str)                 
                newfiles=os.getcwd() + '\\'+ str(i)
                newfilelis.append(i)

        for i in newfilelis:
            nedi=i[:-4]+ '\\'+ str(i)
            cmd_01='copy /Y %s %s'%(str(i),str(nedi))
            print(cmd_01)
            res=os.popen(cmd_01)
            output_str=res.read()
            print(output_str)

        direlist=[dirs for root,dirs,files in os.walk(os.getcwd())]
        direlist=[ i for i in direlist if i][0]
        print(direlist) 
        for i in direlist:
            os.chdir(os.getcwd() + '\\'+ str(i))              

            Ncores=64
            Number_of_cores=Ncores
            Number_of_nodes=Ncores//64

            inputline=[]
            line1_1='#!/bin/bash'
            inputline.append(line1_1)
            line1_2='#SBATCH -p amd_256 '
            inputline.append(line1_2)  
            line2_1='#SBATCH -N %s'%str(Number_of_nodes)
            inputline.append(line2_1)
            line2_1='#SBATCH -n %s'%str(Number_of_cores)
            inputline.append(line2_1)
            line2_2='#SBATCH -J wcl'
            inputline.append(line2_2)
            line2_3='source /public3/soft/modules/module.sh'
            inputline.append(line2_3)
            line2_4='module load cp2k/8.1-public3'
            inputline.append(line2_4)
            line2_5='srun -n %s cp2k.popt -i %s.inp'%(str(Number_of_cores),i)
            inputline.append(line2_5)

            print('                 The number of cores is %s'%str(Number_of_cores))

            ##############################################################################################################
    #WRITE
            filename = str('%s.sh'%str(i))

            s='\n'.join(inputline)
            fp=open(filename,'wb')
            fp.write(s.encode('utf-8'))
            fp.close()
            os.chdir(self.dire)
       

######################################################################################################################
#批处理脚本 
        inputline=[]
        line1_1='#!/bin/bash'
        inputline.append(line1_1)
        line1_2='##This script can submit all tasks under the same level directory '
        inputline.append(line1_2)
        line2_1='cd .'
        inputline.append(line2_1)  
        line2_1='for i in `ls -F |grep /`'
        inputline.append(line2_1)
        line2_1='do'
        inputline.append(line2_1)
        line2_2='cd $i'
        inputline.append(line2_2)
        line2_3='sbatch *.sh'
        inputline.append(line2_3)
        line2_4='cd ..'
        inputline.append(line2_4)
        line2_5='done'
        inputline.append(line2_5)

######################################################################################################################
#WRITE
        filename = str('batch.sh')

        s='\n'.join(inputline)
        fp=open(filename,'wb')
        fp.write(s.encode('utf-8'))
        fp.close()


    def phono_band(self,filename='aimd.cif'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import phonopy
        import subprocess
        from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
        os.chdir(self.dire)

######################################################################################################################
#####################################################################################################################
# perform phonopy

        direlist=[dirs for root,dirs,files in os.walk(os.getcwd())]
        direlist=[ i for i in direlist if i][0]
        # print(direlist) 
        for i in direlist:
            os.chdir(os.getcwd() + '\\'+ str(i))
            direlist=[files for root,dirs,files in os.walk(os.getcwd())]
            direlist=[ i for i in direlist if i][0]
            match_0=[ i for i in direlist if i.endswith(".xyz")]
            if not match_0:
                print('\n\n ************\n%s is blank in terms of force\n************ \n\n'%str(i))
            match=[ i for i in direlist if i.endswith(".xyz")][0]


                
            nedi='..'+ '\\'+ match
            cmd_01='copy /Y %s %s'%(str(match),str(nedi))

            res=os.popen(cmd_01)
            output_str=res.read()
            print(output_str)
            os.chdir(self.dire)
######################################################################################################################
# perform phonopy
        direlist=[files for root,dirs,files in os.walk(os.getcwd())]
        direlist=[ i for i in direlist if i][0]
        match=[ i for i in direlist if i.endswith(".xyz")]
        cmd=['phonopy','--cp2k','-f' ] 
        for i in match: 
            cmd.append(i)
        res=subprocess.run(cmd)
        print(res)  

        path = [[[0, 0, 0], [0.5, 0, 0.5], [0.625, 0.25, 0.625],[0.375, 0.375, 0.75], [0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.25, 0.75]]]
        labels = [r"$\Gamma$", "X", "U", "K", r"$\Gamma$", "L", "W"]
        qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=80)
        phonon = phonopy.load( "phonopy_disp.yaml",'FORCE_SETS',)
        
        phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
        path_connections=connections
        band_dict=phonon.get_band_structure_dict()

        frequencies=band_dict['frequencies']
        distances=band_dict['distances']
        qpoints=band_dict['qpoints']
        def set_xscale_from_data(frequencies, distances):
            """Set xscale from data."""
            max_freq = max([np.max(fq) for fq in frequencies])
            max_dist = distances[-1][-1]
            xscale = max_freq / max_dist * 1.5
            return xscale



        xscale=set_xscale_from_data(frequencies,distances)
        distances_scaled = [d * xscale for d in distances]

        fig,ax = plt.subplots(figsize=(17.5*self.inch_cm,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)
        count = 0
        for i, (d, f, c) in enumerate(
            zip(distances_scaled, frequencies, connections)
        ):
            if i == 0 and labels is not None:
                curves = ax.plot(d, f, c='grey',linestyle='-',linewidth=1.5,alpha=0.85)
                # curves = ax.plot(d, f, c='grey',linestyle='-',linewidth=5.5,alpha=0.35)
            else:
                ax.plot(d, f,c='grey',linestyle='-',linewidth=1.5,alpha=0.85)
                # ax.plot(d, f,c='grey',linestyle='-',linewidth=5.5,alpha=0.35)
            if not c:
                count += 1        

        lefts = [0]
        rights = []
        for i, c in enumerate(path_connections):
            if not c:
                lefts.append(i + 1)
                rights.append(i)
        seg_indices = [list(range(lft, rgt + 1)) for lft, rgt in zip(lefts, rights)]
        special_points = []
        for indices in seg_indices:
            pts = [distances_scaled[i][0] for i in indices]
            pts.append(distances_scaled[indices[-1]][-1])
            special_points.append(pts)

        l_count = 0
        for spts in special_points:
            ax.xaxis.set_ticks_position("both")
            ax.yaxis.set_ticks_position("both")
            ax.xaxis.set_tick_params(which="both", direction="in")
            ax.yaxis.set_tick_params(which="both", direction="in")
            ax.set_xlim(spts[0], spts[-1])
            ax.set_xticks(spts)
            if labels is None:
                ax.set_xticklabels(
                    [
                        "",
                    ]
                    * len(spts)
                )
            else:
                ax.set_xticklabels(labels,fontsize=8,family='Arial')
                l_count += len(spts)
            ax.plot([spts[0], spts[-1]], [0, 0], linestyle="-.", linewidth=0.8, color="b")

        ax.set_ylabel(r"Frequency (THz)",fontsize=8,family='Arial')
        ax.tick_params(labelsize=8)

        plt.savefig('band.tiff', bbox_inches='tight',transparent=True,format='tiff')#指定分辨率,边界紧，背景透明
        plt.show()
        print ('congratulate！！！')        

    def phono_thermol(self,filename='aimd.cif'):
        import os
        import matplotlib.pyplot as plt
        import pandas as pd
        from scipy.optimize import curve_fit
        import numpy as np
        from scipy.signal import savgol_filter
        from scipy.interpolate import make_interp_spline
        from scipy.interpolate import interp2d
        from scipy import interpolate
        from BaselineRemoval import BaselineRemoval     
        from sklearn.preprocessing import MinMaxScaler
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Lattice, Structure, Molecule
        # from Bio.PDB.PDBParser import PDBParser
        # from Bio.PDB.MMCIFParser import MMCIFParser
        # from Bio.PDB.PDBIO import PDBIO
        import re
        import phonopy
        from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
        import subprocess

        os.chdir(self.dire)

#####################################################################################################################
# perform phonopy

        direlist=[dirs for root,dirs,files in os.walk(os.getcwd())]
        direlist=[ i for i in direlist if i][0]
        # print(direlist) 
        for i in direlist:
            os.chdir(os.getcwd() + '\\'+ str(i))
            direlist=[files for root,dirs,files in os.walk(os.getcwd())]
            direlist=[ i for i in direlist if i][0]

            match=[ i for i in direlist if i.endswith(".xyz")][0]

                
            nedi='..'+ '\\'+ match
            cmd_01='copy /Y %s %s'%(str(match),str(nedi))

            res=os.popen(cmd_01)
            output_str=res.read()
            print(output_str)
            os.chdir(self.dire)
######################################################################################################################
# perform phonopy
        direlist=[files for root,dirs,files in os.walk(os.getcwd())]
        direlist=[ i for i in direlist if i][0]
        match=[ i for i in direlist if i.endswith(".xyz")]
        # print(match)
        cmd=['phonopy','--cp2k','-f' ] 
        for i in match: 
            cmd.append(i)
        res=subprocess.run(cmd)
        print(res) 


        phonon = phonopy.load( "phonopy_disp.yaml",'FORCE_SETS',)


        phonon.run_mesh([20, 20, 20])
        phonon.run_thermal_properties(t_step=10,
                                      t_max=1000,
                                      t_min=0)
        tp_dict = phonon.get_thermal_properties_dict()
        temperatures = tp_dict['temperatures']
        free_energy = tp_dict['free_energy']
        entropy = tp_dict['entropy']
        heat_capacity = tp_dict['heat_capacity']
        fig,ax = plt.subplots(figsize=(17.5*self.inch_cm,(17.5*self.inch_cm)*0.5*0.75),constrained_layout=True,dpi=600)

######################################################################################################################
#WRITE  
        inputline=[]
        for a,b,c,d in zip(temperatures,free_energy,entropy,heat_capacity):
            line=str(a)+'   '+str(b)+'   '+str(c)+'   '+str(d)
            inputline.append(line)
        filename = str('thermol')
        s='\n'.join(inputline)
        fp=open(filename,'w+')
        fp.write(s)
        fp.close()
        
        # plt.axvline(x=930,  c='red',ls='dashdot',lw=1)
        # plt.axvline(x=870,  c='red',ls='dashdot',lw=1)

  
        indexs=[0,1,2,3]
        labels=['Free energy (kJ/mol)', 'Entropy (J/K/mol)', 'Heat capacity (J/K/mol)']
        colors=['red', 'blue', 'black']
        lsfs=list(tp_dict.values())
        for t, F, S, cv in zip(temperatures, free_energy, entropy, heat_capacity):
            print(("%12.3f " + "%15.7f" * 3) % ( t, F, S, cv ))
        for i ,j in enumerate(lsfs[1:]):
            ax.plot(lsfs[0],j, c=colors[i],linestyle='-',linewidth=1.5, label=labels[i])
        font1 = {'family' : 'Arial',
        'weight' : 'normal','size':'8'
        }
        # ax.set_ylim((-3,80))
        # ax.set_xlim((4000,300))
        ax.legend(prop=font1,frameon=False,ncol=1)
        ax.tick_params(labelsize=8)        
        # ax.set_xlabel(r'Wavenumber (cm $^{-1}$)',fontsize=8,family='Arial')
        ax.set_xlabel('Temperatures (K)',fontsize=8,family='Arial')
        # ax.legend(edgecolor='none',facecolor='none', prop=font1,loc='best')

        # ax.set_yticks([])
        # ax.set_xticks([])
        ax.tick_params(labelsize=8)

        plt.savefig('thermol.tiff', bbox_inches='tight',transparent=True,format='tiff')#指定分辨率,边界紧，背景透明
        plt.show()
        print ('congratulate！！！')

