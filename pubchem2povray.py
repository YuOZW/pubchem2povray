import os
import random
import shutil
import numpy as np
from scipy.spatial.transform import Rotation
import scipy.optimize as optimize
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import subprocess
import parameter as param
from PIL import Image, ImageFilter, ImageOps

NUM_OUTPUT_IMAGE = 4

def get_compound_data():
	# try searching compounds 100 times
	for _ in range(100):
		target_mol_id = random.randint(1, 110351373)
		c = pcp.get_compounds(target_mol_id, 'cid')
		if len(c) == 0: continue
		c0 = c[0]
		c0_id = c0.cid
		c0_mf = c0.molecular_formula
		c0_smiles = c0.canonical_smiles
		if c0_id is None: continue
		if c0_mf is None: continue
		if c0_smiles is None: continue
		if c0.covalent_unit_count != 1: continue
		break

	if len(c) == 0:
		return None
	else:
		return c0_id, c0_mf, c0_smiles

def calc_scale(coords):
	coords = coords - np.average(coords, axis=0)
	xmax, xmin = coords[:,0].max(), coords[:,0].min()
	xmax = max(np.abs(xmax), np.abs(xmin))
	ymax, ymin = coords[:,1].max(), coords[:,1].min()
	ymax = max(np.abs(ymax), np.abs(ymin))
	scale = max(xmax/16, ymax/9)
	return scale

def rotate_coord(coords):
	def calc_score(para):
		rot = Rotation.from_rotvec(para)
		coords_new = coords.copy()
		for i in range(len(coords)):
			coords_new[i] = rot.apply(coords[i])

		dist_mat = np.zeros((len(coords), len(coords)), dtype=float)
		for i in range(len(coords)):
			for j in range(i):
				dist_mat[i][j] = np.linalg.norm(coords_new[i][:2] - coords_new[j][:2])

		return 1/np.sum(dist_mat * dist_mat)

	rmse = optimize.basinhopping(calc_score, [0.0, 0.0, 0.0], minimizer_kwargs={"method": "BFGS"})

	rot = Rotation.from_rotvec(rmse.x)
	coords_new = coords.copy()
	for i in range(len(coords)):
		coords_new[i] = rot.apply(coords[i])

	return coords_new


def write_povray(natom, labels, coords, scale, bg_info):
	shutil.copy('data/header.pov', 'temp.pov')
	with open('temp.pov', 'a') as fout:
		fout.write(f'{bg_info}\n\n')
		fout.write('merge{\n')
		for i in range(natom):
			x, y, z = coords[i][0]/10, coords[i][1]/10, coords[i][2]/10
			r, g, b = param.color[labels[i]][0], param.color[labels[i]][1], param.color[labels[i]][2]
			line = 'sphere {<' + f'{x:.5f},{y:.5f},{z:.5f}' + '>, '
			if labels[i] == 'H':
				line = line + '0.0252 '
			else:
				line = line + '0.0429 '
			line = line + 'pigment {rgbt<' + f'{r:.3f},{g:.3f},{b:.3f},0.000' + '>} }\n'
			fout.write(line)
		fout.write('\n')

		for i in range(natom):
			for j in range(i):
				length = np.linalg.norm(coords[i]-coords[j])
				length_threshold = param.param[labels[i]][0] + param.param[labels[j]][0]
				if length < length_threshold * 1.5:
					x, y, z = coords[i][0]/10, coords[i][1]/10, coords[i][2]/10
					line = 'cylinder {<' + f'{x:.5f},{y:.5f},{z:.5f}' + '>, '
					x, y, z = coords[j][0]/10, coords[j][1]/10, coords[j][2]/10
					line = line + '<' + f'{x:.5f},{y:.5f},{z:.5f}' + '>, 0.0076 open pigment{rgbt<0.800,0.800,0.800,0.000>} }\n'
					fout.write(line)
		fout.write(f'scale {0.4/scale:.5f}\n')
		fout.write('no_shadow\n')
		fout.write('}\n')
		fout.write('\n')

		fout.write('blob{\n')
		for i in range(natom):
			x, y, z = coords[i][0]/10, coords[i][1]/10, coords[i][2]/10
			r, g, b = param.color[labels[i]][0], param.color[labels[i]][1], param.color[labels[i]][2]
			line = 'sphere {<' + f'{x:.5f},{y:.5f},{z:.5f}' + '>, '
			line = line + str(param.param[labels[i]][1] / 8) + ', 1 '
			line = line + 'pigment {rgbt<' + f'{r:.3f},{g:.3f},{b:.3f},0.500' + '>} }\n'
			fout.write(line)
		fout.write('threshold 0.5\n')
		fout.write(f'scale {0.4/scale:.5f}\n')
		fout.write('no_shadow\n')
		fout.write('}\n')
		fout.write('\n')


def postprocess(image_name):
	if not os.path.isfile(f'{image_name}_1.png'):
		return False
	im1 = np.array(Image.open(f'{image_name}_1.png'))

	if not os.path.isfile(f'{image_name}_2.png'):
		return False
	im2 = np.array(Image.open(f'{image_name}_2.png'))

	im = im1.copy()
	im[:,:,0] = im2[:,:,0]
	im[:,:,1] = im1[:,:,1]
	im[:,:,2] = im1[:,:,2]
	im = Image.fromarray(im)

	mask1 = np.abs(im1[:,:,0] - im2[:,:,0]).astype(np.uint8)
	mask2 = np.abs(im2[:,:,1] - im1[:,:,1]).astype(np.uint8)
	mask = np.stack([mask1, mask2], 2)
	mask = np.amin(mask, axis=2).astype(np.uint8)

	mask = Image.fromarray(mask)
	mask = mask.filter(ImageFilter.MedianFilter())
	mask = ImageOps.invert(mask)
	im.putalpha(mask)

	im.save(f'{image_name}.png')
	if not os.path.isfile(f'{image_name}.png'):
		return False

	return True


def generate_mol_image(image_name):
	# get compound data
	c0_id, c0_mf, c0_smiles = get_compound_data()
	if len(c0_mf) > 50:
		c0_mf = c0_mf[:47] + '...'
	compound_description = f'[CID: {c0_id}] {c0_mf}'

	# first confomation search using MMFF
	m = Chem.MolFromSmiles(c0_smiles)
	m = Chem.AddHs(m)
	cids = AllChem.EmbedMultipleConfs(m, numConfs=100, pruneRmsThresh=2, numThreads=0)
	prop = AllChem.MMFFGetMoleculeProperties(m)
	energies = []
	for i in cids:
		ff = AllChem.MMFFGetMoleculeForceField(m, prop, confId=i)
		if ff is None:
			return False
		ff.Minimize()
		energies.append((ff.CalcEnergy(), i))
	energies.sort()

	# write stable conformer
	writer = Chem.SDWriter("temp.sdf")
	writer.write(m, confId=energies[0][1])
	writer.close() 

	# write charge and spin state for second confomation search using GFN2-xTB
	charge = Chem.GetFormalCharge(m)
	with open('.CHRG', 'w') as fout: fout.write(f'{charge}')
	uhf = Descriptors.NumRadicalElectrons(m)
	with open('.UHF', 'w') as fout: fout.write(f'{uhf}')

	# run xTB
	if os.path.isfile('xtbopt.sdf'):
		os.remove('xtbopt.sdf')
	subprocess.run(['sh', 'run_xtb.sh'])
	if not os.path.isfile('xtbopt.sdf'):
		return False

	# read structure
	labels = []
	coords = []
	bonds = []
	with open('xtbopt.sdf', 'r') as fin:
		fin.readline(); fin.readline(); fin.readline()
		line = fin.readline()
		natom = int(line[0:3])
		nbond = int(line[3:6])
		for i in range(natom):
			line = fin.readline()
			labels.append(line[30:33].strip())
			coords.append([float(line[0:10]), float(line[10:20]), float(line[20:30])])
		for i in range(nbond):
			line = fin.readline()
			bonds.append([int(line[0:3])-1, int(line[3:6])-1])

	# calculate xmax and ymax
	coords = np.array(coords)
	coords = rotate_coord(coords)
	scale = calc_scale(coords)

	# run Pov-ray
	if os.path.isfile('temp.png'):
		os.remove('temp.png')

	write_povray(natom, labels, coords, scale, 'background{ rgb<1.000,0.000,0.000> }')
	subprocess.run(['sh', 'run_povray.sh'])
	if not os.path.isfile('temp.png'):
		return False
	subprocess.run(['mv', 'temp.png', f'{image_name}_1.png'])

	write_povray(natom, labels, coords, scale, 'background{ rgb<0.000,1.000,0.000> }')
	subprocess.run(['sh', 'run_povray.sh'])
	if not os.path.isfile('temp.png'):
		return False
	subprocess.run(['mv', 'temp.png', f'{image_name}_2.png'])

	# finalize
	with open(f'{image_name}.txt', 'w') as fout:
		fout.write(compound_description)

	res = postprocess(image_name)
	if res == False:
		return False

	return True


if __name__ == '__main__':
	imege_id = 0
	for _ in range(100):
		res = generate_mol_image(f'image{imege_id+1}')
		if res == True:
			print(f'[{imege_id+1}/{NUM_OUTPUT_IMAGE}] Done')
			imege_id += 1
		if imege_id >= NUM_OUTPUT_IMAGE:
			break