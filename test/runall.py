# -*- coding: utf-8 -*-
import ROOT as r
from multiprocessing import Pool
import sys
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".x loadTPGSimAnalysis.C")

if len(sys.argv) == 2: nCores = int(sys.argv[1])
else:                  nCores = 1

#suffix = ""
suffix = "_AllCuts"
#suffix = "_newDPhiDEtacuts"
#suffix = "_newSL2requirement"
#suffix = "_newTimecut"

#path = "/pool/ciencias/userstorage/sscruz/DT_TDR_12sep/%s.root"
#path = "/pool/ciencias/userstorage/vrbouza/ntuples/DT_ntuples/shifts_offline_2019/2019_11_05/%s.root"
path = "/pool/ciencias_users/user/vrbouza/www/Miscelánea/2019_11_06_temporalntuples/%s.root"


#def ExecuteTheThing(pu, age, rpc, qual, ind = -99, indmax = +99):
def ExecuteTheThing(tsk):
    pu, age, rpc, qual, ind, indmax, nu = tsk
    inname  = '{neutrinos}{pileup}_{agemaster}_{rpcs}{scenario}'.format(neutrinos = "nu_" * nu + "",
                                                                        pileup = pu,
                                                                        agemaster = age if "noage" in age else "age",
                                                                        rpcs = rpc,
                                                                        scenario = ("_youngseg_" + age) if "noage" not in age else "")
    outname = ("results_eff{neutrinos}_{pileup}_{agemaster}_{rpcs}_{quality}{scenario}{sfx}.root".format(neutrinos = "nu_" * nu + "",
                                                                                                    pileup = pu,
                                                                                                    agemaster = age if "noage" in age else "age",
                                                                                                    rpcs = rpc,
                                                                                                    quality = qual,
                                                                                                    scenario = ("_" + age) if "noage" not in age else "",
                                                                                                    sfx = suffix)
               * (ind == -99 and indmax == +99)
               + "results_eff{neutrinos}_{pileup}_{agemaster}_{rpcs}_{quality}_{indexes}{scenario}{sfx}.root".format(neutrinos = "_nu" * nu + "",
                                                                                                                pileup = pu,
                                                                                                                agemaster = age if "noage" in age else "age",
                                                                                                                rpcs = rpc,
                                                                                                                quality = qual,
                                                                                                                indexes = "ind%i"%(ind) * (ind != -99) + "indmax%i"%(indmax) * (indmax != +99),
                                                                                                                scenario = ("_" + age) if "noage" not in age else "",
                                                                                                                sfx = suffix)
               * (ind != -99 or indmax != +99))

    print "> Analysing file", inname
    analysis = r.DTNtupleTPGSimAnalyzer(path%inname, "./results/" + outname, qual, ind, indmax)
    analysis.Loop()
    print "> File saved as", outname
    return


tasks = []
for neutr in [False]:
    for pu in ["nopu", "pu140", "pu200", "pu250"]:
        for age in ["noage", "muonage_norpcage_nofail_3000_OLD", "muonage_norpcage_fail_3000", "muonage_norpcage_nofail_3000", "muonage_norpcage_fail_1000", "muonage_norpcage_nofail_1000"]:
            for rpc in ["withrpc", "norpc"]:
                #for qual in ["", "nothreehits", "higherthanfour", "higherthanfourvetoing", "qualityOR"]:
                for qual in ["qualityOR"]:
                #for qual in 'higherthanfour,higherthanfourvetoing'.split(','):
                    if (not neutr) and (("140" in pu) or ("250" in pu) or ("300" in pu)): continue
                    tasks.append( (pu, age, rpc, qual, -99, +99, neutr) )
                    if qual == "":
                        tasks.append( (pu, age, rpc, qual, 0, +99, neutr) )
                        tasks.append( (pu, age, rpc, qual, -99, 2, neutr) )

#tasks.append( ("nopu", "age", "withrpc", "higherthanfourvetoing", -99, +99) )


if nCores == 1:
    for task in tasks: ExecuteTheThing(task)
else:
    print "\n> Launching plotting processes..."
    pool = Pool(nCores)
    pool.map(ExecuteTheThing, tasks)
    pool.close()
    pool.join()
    del pool


