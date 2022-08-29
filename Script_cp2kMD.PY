import sys
sys.path.append(r"E:\Desktop\VASP practical\Pymatgen")
from cp2kMD import cp2kMD

my_cp2kMD = cp2kMD(r'F:\cp2k practical\SFeS\110_dio')


'''
结果检查
'''
# my_cp2kMD.POSTAIMD(filename='run-NH2O-MD-T300.xyz',rdf_fir='type O',rdf_sec='type H')

# my_cp2kMD.Energychek()

# my_cp2kMD.forcechek(atomindex=15)

# my_cp2kMD.Ensemblechek('slurm-3354687.out1')

'''
热力学量计算
'''
# my_cp2kMD.freeenergy()

# my_cp2kMD.IRspectrum(corr=300,FWHM=60,xlim_1=1660,xlim_2=1200)#FWHM=34


# sulf_ener={(0, 0, 2):2.8, 
#                         (2, 0, 0): 4.8}
# my_cp2kMD.Wullff(surface_energy=sulf_ener,directions=(3,4,4),x=0.47,y=-16)


# my_cp2kMD.phono_pre()
# my_cp2kMD.phono_band()
# my_cp2kMD.phono_thermol()

'''
波函数分析
'''
# my_cp2kMD.pdos(nelement=1)#pdos,nelement为元素顺序

# my_cp2kMD.lpdos(nlist=1)

# my_cp2kMD.moldenformultiwfn()#将cp2k生成的molden文件转化为multiwfn可读的格式

# my_cp2kMD.IGA()#本机计算很慢，这是生成用LINUX算IGA的相关文件,需放入上一个函数产生的molden文件

# my_cp2kMD.den_dif(vdw_atom=4,atomindex=83)#算两个特定原子附近的差分cube，请将被减数名为minuend.cube，减数名为subtraction.cube

# my_cp2kMD.vmdrender()#将VMD生成的vmdscene.dat渲染到指定目录下

# my_cp2kMD.localizecube(filename='func1.cub',vdw_atom=1.5,atomindex=118)#atomindex=118

# my_cp2kMD.Elec_poten_peri()

# my_cp2kMD.Elec_potential()

'''
输入文件生成
'''
my_cp2kMD = cp2kMD(r'F:\cp2k practical\SFeS\SFeS_S4_01')
# # # # # # # #############提供的cif文件名为aimd####################################
# # # # # # # ##################################################################################
# my_cp2kMD.pdbtocif()

my_cp2kMD.frame(sccs=0,semethod='QS',runtype='MD')#FIST,QS,QMMM,MD,CELL_OPT,VIBRATIONAL_ANALYSIS,ENERGY #

my_cp2kMD.motiongenerate(ensemble='NPT_F',fixatoms=(0,10000))#NVT,NVE,NPT_F#fixatoms=['ssS','Oss','Cssa']/(1.5,23)z坐标在两个数之间是自由的

my_cp2kMD.dftgenerate(method='DI',ldoslist=[1,2])#method='DI'/'OT'[179,88,98,21,23]

my_cp2kMD.colvargenerate(atomlist=[[1,2],[1,2]])

my_cp2kMD.pdbgenerate()

my_cp2kMD.shgenerate(Ncores=64)


# my_cp2kMD.mmgenerate()