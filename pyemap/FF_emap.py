from .aminoacids import params_dict 
import numpy as np
from .data import res_name_to_char, side_chain_atoms, char_to_res_name, site_energies
from .SASA_data import cation_sasa
from .SASA_data import anion_sasa
from Bio.PDB.DSSP import DSSP
import warnings
import traceback
import re







def create_forcefield(inputpdb,emap):
	natoms=0
	with open(inputpdb, 'r') as pdb:
		for line in pdb:
			resname=line.split()[3]
			if (line.split()[0] == 'ATOM'):
				natoms+=1


	reslist = []
	atomlist=[]
	moleculelist=[]

	ppdb = open(inputpdb, 'r')
	data = open(inputpdb).read()
	with open(inputpdb, 'r') as pdb:
		for line in pdb:
			if (line.split()[0] == 'ATOM'):



			#	print('all atoms', line.split()[1],line.split()[2],line.split()[3],line.split()[5],line.split()[6],line.split()[7],line.split()[8])

			#	print(len(params_dict[line.split()[3]]))
			
				checkstring=f"{line.split()[3]:<3}{line.split()[4]:^3}{line.split()[5]:>3}"
				if checkstring in line:
					print(line)
				print('building file', data.count(checkstring))
				print(checkstring)
				print((len(params_dict[line.split()[3]])))
				print((params_dict[line.split()[3]]))
			#	print('checkstring',checkstring)
			#	print('checking for length', len(params_dict[line.split()[3]]), data.count(checkstring))
				
				if (len(params_dict[line.split()[3]]) != data.count(checkstring)):
					#print('continued', params_dict[line.split()[3]])
					b=0
				else:
					#print(len(params_dict[line.split()[3]]))

				
					ff_param = params_dict[line.split()[3]][line.split()[2]]
					#print('charge',line.split()[1],line.split()[2],line.split()[3], line.split()[5], ff_param)
				#print(ff_param)
				#print(line.split()[2],ff_param)
					atom = (line.split()[6], line.split()[7], line.split()[8], line.split()[2], line.split()[3], ff_param, line.split()[4], line.split()[5], line.split()[1] )
				
					moleculelist.append(atom) 

					atomlist.append(line.split()[1])



					if line.split()[5] not in reslist: 
			#	if line.split()[5] =='13' && line.split()[5] not in reslist:
						reslist.append(line.split()[5])



	Fields=[0]*len(reslist)*2
	#print(len(reslist))

	surface_numbers=[]


	for i in range(0,len(emap.surface_residues)):
			surface_num=re.split('(\d+)',emap.surface_residues[i])[1]
			surface_numbers.append(surface_num)
	print('surface numbers', surface_numbers)

	for residue_number in reslist:
		if residue_number in surface_numbers:
			Fields[int(residue_number)]=0.00
		else:
		

			Field = build_E_field(residue_number,moleculelist,natoms,resname)
	#	print('residue_number', residue_number, field)
		#print('index', residue_number, type(residue_number))
			residue_number = int(residue_number)
		#print('residue number', residue_number)
		#print(Fields)
			Fields[residue_number] = Field

	return(Fields)



def build_SASA(filename,model):

    cation_residues=['V', 'A', 'G', 'L', 'I', 'M', 'P', 'S','T', 'C', 'N', 'Q', 'D', 'E', 'R', 'K']
    anion_residues=['W', 'F', 'H', 'Y' ]

    """
    Only standard protein residues are currently supported. Non-protein and user specified custom residues cannot be
    classified as surface exposed using this criteria.

    Parameters
    ---------
    filename: str
        Name of pdb file to be analyzed
    model: :class:`Bio.PDB.Model.Model`
        Model under analysis
    node_list : list of str
        List containing which standard residues are included in analysis

    References
    ---------
    Tien, M. Z.; Meyer, A. G.; Sydykova, D. K.; Spielman, S. J.; Wilke, C. O. PLoS ONE 2013, 8 (11).
        Reference for relative solvent accessibility cutoff of 0.05, and for MaxASA values
    """
    solv_cont ={'cation':{},'anion':{}}
    try:
        dssp = DSSP(model, filename, acc_array="Wilke")
        for key in dssp.keys():

            goal_str = dssp[key][1] + str(key[1][1]) + "(" + str(key[0]) + ")"
           # print('key1', key[1][1], 'dssp', dssp[key][3], 'goal_str', goal_str)
           # print(dssp[key])
           # print(dssp[key][1]) #C


            #this part has nothing to do with dssp and is a printing checkstep for me
      
            if dssp[key][1] in cation_sasa.keys():

            		SASA_hole = (1-float(dssp[key][3])) * (float(cation_sasa[dssp[key][1]]))
            		solv_cont['cation'][goal_str] = SASA_hole
            	#	print('SASA hole',goal_str, SASA_hole)


            if dssp[key][1] in anion_sasa.keys():

            		SASA_electron = float(dssp[key][3]) * (1-(anion_sasa[dssp[key][1]]))
            		solv_cont['anion'][goal_str] = SASA_electron
            	#	print('SASA electron', goal_str, SASA_hole)





    except Exception as e:
    		print('Error',e)
    		print(traceback.format_exc())
    		warnings.warn("Unable to calculate solvent accessibility. Check that DSSP is installed.",RuntimeWarning,stacklevel=2)
   	

    return solv_cont




#


def build_E_field(residue_number,moleculelist,natoms,resname):




		charge_total=0.00
		#c= 10**-10/(4* np.pi * 8.85418782 *10**-12 )
		c=1
		#print(c)
		com_x, com_y, com_z, res_name=com_res(residue_number,moleculelist)
		adjusted_xfield_sum=0
		adjusted_yfield_sum=0
		adjusted_zfield_sum=0
		xfield_sum=0
		yfield_sum=0
		zfield_sum=0
		shift=0.0
		current_res_shift=[0]*500




		for atom in moleculelist:
			charge_total= float(atom[5])+charge_total
			#print('running charge',float(atom[5]),charge_total)
		#print('length',len(moleculelist))
	#	print('charge total pre', charge_total)
		#print('len mol list', len(moleculelist))
		for atom in moleculelist:
			xfield=0
			yfield=0
			zfield=0
			total_charge=0.0
			adjusted_zfield=0
			adjusted_xfield =0
			adjusted_yfield=0


	
			

			if atom[7] != residue_number:  
			#	print('check', atom[7], residue_number)


				rx = (com_x-float(atom[0]))*1.88973  #I think it is here. boo 
				ry= (com_y-float(atom[1]))*1.88973
				rz = (com_z-float(atom[2]))*1.88973

	
				charge=(float(atom[5]))
		#		print('charge', float(atom[5]), 'charge total', charge_total/len(moleculelist))
				adjusted_charge = float(atom[5])  - charge_total/(len(moleculelist))
			#	print(charge, charge_total/len(moleculelist))

				total_charge += adjusted_charge
				#print('running adjusted charge', total_charge)


				#print(int(residue_number))
				current_res_shift[int(atom[7])] += adjusted_charge/np.linalg.norm([rx,ry,rz])
				


				shift += adjusted_charge/np.linalg.norm([rx,ry,rz])
			#	print('shift', shift)

		#for i in range (1,int(atom[7])):

		#	print(residue_number,27.211382543519*current_res_shift[i], i)

	#	print('total charge', total_charge)
		#print('Current res shift',current_res_shift)
		print('Total shift', residue_number,resname, 27.211382543519*shift)
		return( 27.211382543519*shift)



def com_res(residue_number,moleculelist):
	mass_sum=0
	com_x =0
	com_y =0
	com_z = 0
	molecule=[]
	sca=[]
	x_wsum =0
	y_wsum =0
	z_wsum =0
	for atom in moleculelist:
	#	print(molecule)
	#	print('creating com')
		if (atom[7] ==residue_number) and (atom[4] in side_chain_atoms):
		#	print(atom[7], residue_number, atom[4])
			sca = side_chain_atoms[atom[4]]
			#print('sca',sca)
			if atom[3] in sca:
			#	print('side chain atoms', atom[3])

				mass=find_mass(atom)
				x = float(atom[0])
				y =float(atom[1])
				z = float(atom[2])
				#print('ind coords',atom[0],atom[1], atom[2])
				#print(atom[0], atom[1], atom[2])
				x_wsum += x * mass
				y_wsum += y * mass
				z_wsum += z * mass
			#	print('com_build',x,y,z,mass,atom[3],atom[4],residue_number)
				mass_sum += mass


	com_x = x_wsum / mass_sum
	com_y = y_wsum / mass_sum
	com_z = z_wsum / mass_sum
	#print('com',residue_number, com_x,com_y,com_z)
#	print('com coords', com_x, com_y, com_z)
	return(com_x, com_y, com_z, atom[4])

def find_mass(atom):
	a = atom[3]
	lista=list(a)
	if lista[0] =='N' :
		mass = 14
	if lista[0] == 'C' :
		mass = 12.01
	if lista[0] == 'H' :
		mass = 1.001
	if lista[0] == 'S':
		mass = 32.96
	if lista[0] == 'O':
		mass = 16.02

	return(mass)




		




#print(dir(chargelist))
#print(chargelist)

         
# both objects have different self which
# contain their attributes
    # same output as car.show(audi)

