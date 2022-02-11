import os
import sys
from collections import namedtuple

ScoreData = namedtuple('ScoreData', 'total_score score dslf_fa13 fa_atr fa_dun fa_elec fa_intra_rep fa_intra_sol_xover4  fa_rep fa_sol hbond_bb_sc hbond_lr_bb hbond_sc hbond_sc_bb linear_chainbreak lk_ball_wtd omega overlap_chainbreak p_aa_pp pro_close rama_prepro ref time yhh_planarity description')


def getListLines(file1):
	listLines = [line.rstrip("\t\n") for line in open(file1, "r")]
	return listLines


def main(inputFolder, outputFile):
	listAllData = []
	for file1 in os.listdir(inputFolder):
		if "6hax_" == file1[0:5]\
		and ".score" == file1[len(file1)-6:len(file1)]:
			file1Path = os.path.join(inputFolder, file1)
			listLines = getListLines(file1Path)
			if 3 != len(listLines):
				break
			line = listLines[2]
			total_score = line.split()[1].strip()
			score = line.split()[2].strip()
			dslf_fa13 = line.split()[3].strip()
			fa_atr = line.split()[4].strip()
			fa_dun = line.split()[5].strip()
			fa_elec = line.split()[6].strip()
			fa_intra_rep = line.split()[7].strip()
			fa_intra_sol_xover4 = line.split()[8].strip()
			fa_rep = line.split()[9].strip()
			fa_sol = line.split()[10].strip()
			hbond_bb_sc = line.split()[11].strip()
			hbond_lr_bb = line.split()[12].strip()
			hbond_sc = line.split()[13].strip()
			hbond_sc_bb = line.split()[14].strip()
			linear_chainbreak = line.split()[15].strip()
			lk_ball_wtd = line.split()[16].strip()
			omega = line.split()[17].strip()
			overlap_chainbreak = line.split()[18].strip()
			p_aa_pp = line.split()[19].strip()
			pro_close = line.split()[20].strip()
			rama_prepro = line.split()[21].strip()
			ref = line.split()[22].strip()
			time = line.split()[23].strip()
			yhh_planarity = line.split()[24].strip()
			description = line.split()[25].strip()
			element = ScoreData( total_score = total_score, score = score, dslf_fa13 = dslf_fa13, fa_atr = fa_atr, fa_dun = fa_dun, fa_elec = fa_elec, fa_intra_rep = fa_intra_rep, fa_intra_sol_xover4 = fa_intra_sol_xover4, fa_rep = fa_rep, fa_sol = fa_sol, hbond_bb_sc = hbond_bb_sc, hbond_lr_bb = hbond_lr_bb, hbond_sc = hbond_sc, hbond_sc_bb = hbond_sc_bb, linear_chainbreak = linear_chainbreak, lk_ball_wtd = lk_ball_wtd, omega = omega, overlap_chainbreak = overlap_chainbreak, p_aa_pp = p_aa_pp, pro_close = pro_close, rama_prepro = rama_prepro, ref = ref, time = time, yhh_planarity = yhh_planarity, description = description )
			listAllData.append(element)
	listAllDataSorted = sorted(listAllData, key = lambda x: x.total_score)
	with open(outputFile, "a") as out:
		line = "total_score,score,dslf_fa13,fa_atr,fa_dun,fa_elec,fa_intra_rep,fa_intra_sol_xover4,fa_rep,fa_sol,hbond_bb_sc,hbond_lr_bb,hbond_sc,hbond_sc_bb,linear_chainbreak,lf_ball_wtd,omega,overlap_chainbreak,p_aa_pp,pro_close,rama_prepro,ref,time,yhh_planarity,description\n"
		out.write(line)
		for el in listAllDataSorted:
			line = el.total_score + "," + el.score + "," + el.dslf_fa13 + "," + el.fa_atr + "," + el.fa_dun + "," + el.fa_elec + "," + el.fa_intra_rep + "," + el.fa_intra_sol_xover4 + "," + el.fa_rep + "," + el.fa_sol + "," + el.hbond_bb_sc + "," + el.hbond_bb_sc + "," + el.hbond_lr_bb + "," + el.hbond_sc + "," + el.hbond_sc_bb + "," + el.linear_chainbreak + "," + el.lk_ball_wtd + "," + el.omega + "," + el.overlap_chainbreak + "," + el.p_aa_pp + "," + el.pro_close + "," + el.rama_prepro + "," + el.ref + "," + el.time + "," + el.yhh_planarity + "," + el.description + "\n"
			out.write(line)




if "__main__" == __name__:
	if 3 != len(sys.argv):
		print("please provide:\n(i). input folder with *.score files;\n(ii). output file to save results;")
		sys.exit()
	inputFolder = sys.argv[1]
	outputFile = sys.argv[2]
	if False == os.path.exists(inputFolder)\
	or False == os.path.isdir(inputFolder):
		print(f'folder {inputFolder} does not exist')
		sys.exit()
	if True == os.path.exists(outputFile):
		os.remove(outputFile)
	main(inputFolder, outputFile)
