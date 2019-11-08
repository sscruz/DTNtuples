# -*- coding: utf-8 -*-
import ROOT as r
from copy import deepcopy
import CMS_lumi, os
r.gROOT.SetBatch(True)

#path = "~sscruz/www/DT_TDR/2019_12_09_plots_eff_withHBaged_noquality/"
#path = "~vrbouza/www/Miscelánea/2019_11_04_plots_eff_shiftsoff/IndividualApplyingOfCuts/"
path = "~vrbouza/www/Miscelánea/2019_11_07_plots_eff_shiftsoff/AgeingComparisons/"
#path = "~vrbouza/www/Miscelánea/2019_11_04_plots_eff_shiftsoff/"

openingpath  = "./results/"
plotscaffold = "hEff_{st}_{al}_{ty}"
savescaffold = "hEff_{pu}{qu}{id}"
chambTag     = ["MB1", "MB2", "MB3", "MB4"]

suffix       = ""
#suffix       = "_AllCuts"
#suffix       = "_newDPhiDEtacuts"
#suffix       = "_newSL2requirement"
#suffix       = "_newTimecut"

#extrasuffix = ""
extrasuffix  = "_AllScenarios"
#extrasuffix  = "_ONLYAMRPCWITHAM"
#extrasuffix  = "_ONLYAM"
#extrasuffix  = "_ONLYHB"
scenariolist = ["noage", "muonage_norpcage_nofail_3000_OLD", "muonage_norpcage_fail_3000", "muonage_norpcage_nofail_3000", "muonage_norpcage_fail_1000", "muonage_norpcage_nofail_1000"]


def makeresplot(hlist, aged, algo, qual = "", pu = "nopu", ind = "", neutr = ""):
    print "\nObtaining intermediate plot for algo", algo, "which has", aged, "ageing and considering", pu

    openfile = "results_eff_" + neutr + pu + "_" + (aged if "noage" in aged else "age") + "_" + ("with" * ("RPC" in algo) + "no" * ("RPC" not in algo)) + "rpc_" + qual + "_" * ("noage" not in aged) + ind + aged * ("noage" not in aged) + suffix + ".root"

    print "The file that it's going to be opened is", openfile
    res = r.TFile.Open(openingpath + openfile, "read")

    print "Then, we're gonna get the histograms labeled as follows:", plotscaffold.format(al = algo.replace("+RPC", ""), st = chambTag[0], ty = "matched")

    hmatched = [res.Get(plotscaffold.format(al = algo.replace("+RPC", ""), st = chambTag[ich], ty = "matched")) for ich in range(4)]
    htotal   = [res.Get(plotscaffold.format(al = algo.replace("+RPC", ""), st = chambTag[ich], ty = "total")) for ich in range(4)]

    savehisto = "hEff_{al}_{ag}_{pileup}".format(al = algo, ag = aged, pileup = pu) + "_{q}".format(q = qual) * (qual != "") + "_{id}".format(id = ind)* (ind != "")

    print "And finally the histogram that is going to be saved will have this name:", savehisto
    resplot = r.TH1D(savehisto, "", 20, -0.5, 19.5)
    
    ibin = 1
    for ich in range(4):
        for iwh in range(1, 6):
            #print "st", ich, "wh", iwh, "valmatched:", hmatched[ich].GetBinContent(iwh), "valtotal:", htotal[ich].GetBinContent(iwh)
            if htotal[ich].GetBinContent(iwh) == 0:
                print "WARNING: denominator for chamber", ich, "and wheel", iwh, "is ZERO: fixing efficiency to zero too."
                effval = 0
            else:
                effval = hmatched[ich].GetBinContent(iwh) / htotal[ich].GetBinContent(iwh)

            if effval > 1: print "WARNING: eff. for chamber", ich, "and wheel", iwh, "is", effval, ". Fixing to 1."

            resplot.SetBinContent(ibin, effval * (effval < 1) + 1 * (effval >= 1))
            eff = r.TEfficiency('kk','',1,-0.5,0.5)
            eff.SetTotalEvents(1, int(htotal[ich].GetBinContent(iwh)))
            eff.SetPassedEvents(1, int(hmatched[ich].GetBinContent(iwh)))
            if (eff.GetEfficiencyErrorLow(1) - eff.GetEfficiencyErrorUp(1)) > 0.05: print 'warning, bin asymmetric'
            resplot.SetBinError( ibin, max(eff.GetEfficiencyErrorLow(1), eff.GetEfficiencyErrorUp(1)))
            del eff
            ibin += 1

    hlist.append(deepcopy(resplot))
    res.Close(); del hmatched, htotal, res, resplot
    return


#lowlimityaxis  = 0.2
lowlimityaxis  = 0.8
#lowlimityaxis  = 0.4
highlimityaxis = 1
markersize     = 1
yaxistitle     = "Efficiency (adim.)"
yaxistitleoffset= 1.5
xaxistitle     = "Wheel"
#legxlow        = 0.3075 + 2 * 0.1975 # Old
#legxlow        = 0.3075 + 1 * 0.1975 # Good
legxlow        = 0.3075
#legxlow        = 0.3075 + 1/3. * 0.1975
legylow        = 0.3 # Good one
#legylow        = 0.25
#legylow        = 0.1 # floor one
legxhigh       = 0.9
legyhigh       = 0.5 # Good one
#legyhigh       = 0.6
#legyhigh       = 0.45 # floor one
legtextsize    = 0.02

markertypedir  = {}
markertypedir["AM_age"]       = 24
markertypedir["AM_noage"]     = 20
#markertypedir["AM+RPC_age"]   = 29
#markertypedir["AM+RPC_noage"] = 29
markertypedir["AM+RPC_age"]   = 26
markertypedir["AM+RPC_noage"] = 22
markertypedir["HB_noage"]     = 21
markertypedir["HB_age"]       = 25

markertypedir["AM_muonage_norpcage_nofail_3000_OLD"]     = markertypedir["AM_age"]
markertypedir["AM+RPC_muonage_norpcage_nofail_3000_OLD"] = markertypedir["AM+RPC_age"]
markertypedir["HB_muonage_norpcage_nofail_3000_OLD"]     = markertypedir["HB_age"]

markertypedir["AM_muonage_norpcage_fail_3000"]           = markertypedir["AM_age"]
markertypedir["AM+RPC_muonage_norpcage_fail_3000"]       = markertypedir["AM+RPC_age"]
markertypedir["HB_muonage_norpcage_fail_3000"]           = markertypedir["HB_age"]

markertypedir["AM_muonage_norpcage_nofail_3000"]         = markertypedir["AM_age"]
markertypedir["AM+RPC_muonage_norpcage_nofail_3000"]     = markertypedir["AM+RPC_age"]
markertypedir["HB_muonage_norpcage_nofail_3000"]         = markertypedir["HB_age"]

markertypedir["AM_muonage_norpcage_fail_1000"]           = markertypedir["AM_age"]
markertypedir["AM+RPC_muonage_norpcage_fail_1000"]       = markertypedir["AM+RPC_age"]
markertypedir["HB_muonage_norpcage_fail_1000"]           = markertypedir["HB_age"]

markertypedir["AM_muonage_norpcage_nofail_1000"]         = markertypedir["AM_age"]
markertypedir["AM+RPC_muonage_norpcage_nofail_1000"]     = markertypedir["AM+RPC_age"]
markertypedir["HB_muonage_norpcage_nofail_1000"]         = markertypedir["HB_age"]


markercolordir  = {}
markercolordir["AM_age"]       = 46
markercolordir["AM+RPC_age"]   = 46
markercolordir["AM+RPC_noage"] = 9
markercolordir["HB_noage"]     = 9
markercolordir["AM_noage"]     = 9
markercolordir["HB_age"]       = 46

markercolordir["AM_muonage_norpcage_nofail_3000_OLD"]     = markercolordir["AM_age"]
markercolordir["AM+RPC_muonage_norpcage_nofail_3000_OLD"] = markercolordir["AM+RPC_age"]
markercolordir["HB_muonage_norpcage_nofail_3000_OLD"]     = markercolordir["HB_age"]

###
#markercolordir["AM_muonage_norpcage_fail_3000"]           = markercolordir["AM_age"]
#markercolordir["AM+RPC_muonage_norpcage_fail_3000"]       = markercolordir["AM+RPC_age"]
#markercolordir["HB_muonage_norpcage_fail_3000"]           = markercolordir["HB_age"]
markercolordir["AM_muonage_norpcage_fail_3000"]           = 4
markercolordir["AM+RPC_muonage_norpcage_fail_3000"]       = 4
markercolordir["HB_muonage_norpcage_fail_3000"]           = 4
###

markercolordir["AM_muonage_norpcage_nofail_3000"]         = 7
markercolordir["AM+RPC_muonage_norpcage_nofail_3000"]     = 7
markercolordir["HB_muonage_norpcage_nofail_3000"]         = 7

markercolordir["AM_muonage_norpcage_fail_1000"]           = 2
markercolordir["AM+RPC_muonage_norpcage_fail_1000"]       = 2
markercolordir["HB_muonage_norpcage_fail_1000"]           = 2

markercolordir["AM_muonage_norpcage_nofail_1000"]         = 6
markercolordir["AM+RPC_muonage_norpcage_nofail_1000"]     = 6
markercolordir["HB_muonage_norpcage_nofail_1000"]         = 6

namedir  = {}
namedir["AM_age"]       = "DT AM w/ ageing"
namedir["AM+RPC_age"]   = "DT AM w/ ageing w/ RPC"
namedir["AM+RPC_noage"] = "DT AM w/ RPC"
namedir["HB_noage"]     = "DT HB"
namedir["AM_noage"]     = "DT AM"
namedir["HB_age"]       = "DT HB w/ ageing"

scenariolist = ["noage", "muonage_norpcage_nofail_3000_OLD", "muonage_norpcage_fail_3000", "muonage_norpcage_nofail_3000", "muonage_norpcage_fail_1000", "muonage_norpcage_nofail_1000"]

namedir["AM_muonage_norpcage_nofail_3000_OLD"]     = "DT AM w/ DT 3000fb-1 OLD"
namedir["AM+RPC_muonage_norpcage_nofail_3000_OLD"] = "DT AM w/ DT 3000fb-1 OLD w/ RPC"
namedir["HB_muonage_norpcage_nofail_3000_OLD"]     = "DT HB w/ DT 3000fb-1 OLD"

### TAN MAL: CAMBIÁOS LOS FAIL Y NO FAIL DE 3000
#namedir["AM_muonage_norpcage_fail_3000"]           = "DT AM w/ DT+RPC 3000fb-1"
#namedir["AM+RPC_muonage_norpcage_fail_3000"]       = "DT AM w/ DT+RPC 3000fb-1 w/ RPC"
#namedir["HB_muonage_norpcage_fail_3000"]           = "DT HB w/ DT+RPC 3000fb-1"

#namedir["AM_muonage_norpcage_fail_3000"]           = namedir["AM_age"]
#namedir["AM+RPC_muonage_norpcage_fail_3000"]       = namedir["AM+RPC_age"]
#namedir["HB_muonage_norpcage_fail_3000"]           = namedir["HB_age"]
#TAN MAL: CAMBIÁOS LOS FAIL Y NO FAIL DE 3000
namedir["AM_muonage_norpcage_fail_3000"]           = "DT AM w/ DT 3000fb-1"
namedir["AM+RPC_muonage_norpcage_fail_3000"]       = "DT AM w/ DT 3000fb-1 w/ RPC"
namedir["HB_muonage_norpcage_fail_3000"]           = "DT HB w/ DT 3000fb-1"
###

#TAN MAL: CAMBIÁOS LOS FAIL Y NO FAIL DE 3000
namedir["AM_muonage_norpcage_nofail_3000"]         = "DT AM w/ DT+RPC 3000fb-1"
namedir["AM+RPC_muonage_norpcage_nofail_3000"]     = "DT AM w/ DT+RPC 3000fb-1 w/ RPC"
namedir["HB_muonage_norpcage_nofail_3000"]         = "DT HB w/ DT+RPC 3000fb-1"

#namedir["AM_muonage_norpcage_nofail_3000"]         = "DT AM w/ DT 3000fb-1"
#namedir["AM+RPC_muonage_norpcage_nofail_3000"]     = "DT AM w/ DT 3000fb-1 w/ RPC"
#namedir["HB_muonage_norpcage_nofail_3000"]         = "DT HB w/ DT 3000fb-1"

namedir["AM_muonage_norpcage_fail_1000"]           = "DT AM w/ DT+RPC 1000fb-1"
namedir["AM+RPC_muonage_norpcage_fail_1000"]       = "DT AM w/ DT+RPC 1000fb-1 w/ RPC"
namedir["HB_muonage_norpcage_fail_1000"]           = "DT HB w/ DT+RPC 1000fb-1"

namedir["AM_muonage_norpcage_nofail_1000"]         = "DT AM w/ DT 1000fb-1"
namedir["AM+RPC_muonage_norpcage_nofail_1000"]     = "DT AM w/ DT 1000fb-1 w/ RPC"
namedir["HB_muonage_norpcage_nofail_1000"]         = "DT HB w/ DT 1000fb-1"



def combineresplots(hlist, qual = "", pued = False, ind = "", zoom = "normal"):
    print "Combining list of plots that has", pued, "pile-up,", qual, "quality and", ind, "index."
    if len(hlist) == 0: raise RuntimeError("Empty list of plots")

    c = r.TCanvas("c", "c", 800, 800)
    c.SetLeftMargin(0.11)
    c.SetGrid()

    leg = r.TLegend(legxlow, legylow, legxhigh, legyhigh)
    #leg.SetTextSize(legtextsize);
    hlist[0].SetStats(False)
    #hlist[0].SetTitle("L1 DT Phase 2 algorithm efficiency comparison")
    hlist[0].GetYaxis().SetRangeUser(lowlimityaxis if (zoom == "normal") else 0.65 if (zoom == "zoomin") else 0.8 if (zoom == "extremezoomin") else 0.4, (highlimityaxis + (highlimityaxis - (lowlimityaxis if (zoom == "normal") else 0.4)) * 0.1) if (zoom != "zoomin" and zoom != "extremezoomin") else 1.0125 if (zoom == "extremezoomin") else 1.05)
    hlist[0].GetYaxis().SetTitleOffset(yaxistitleoffset)
    hlist[0].GetYaxis().SetTitle(yaxistitle)
    hlist[0].GetXaxis().SetTitle(xaxistitle)
    hlist[0].GetXaxis().SetNdivisions(120)
    ilabel = 1
    for ich in range(4):
        for iwh in range(-2, 3):
            hlist[0].GetXaxis().ChangeLabel(ilabel, -1, -1, -1, -1, -1, (iwh > 0) * "+" + str(iwh))
            ilabel += 1

    for iplot in range(len(hlist)):
        hlist[iplot].SetMarkerSize(markersize)    # .replace("_nopu","").replace("_pu","").replace("hEff_","")
        #hlist[iplot].SetMarkerStyle(markertypedir[( hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2] )])
        #hlist[iplot].SetMarkerColor(markercolordir[hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2]])
        ##leg.AddEntry(hlist[iplot], namedir[hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2]], "P")
        #leg.AddEntry(hlist[iplot], namedir[hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2]])

        hlist[iplot].SetMarkerStyle(markertypedir[( hlist[iplot].GetName().split("_pu")[0].split("hEff_")[1] if pued else hlist[iplot].GetName().split("_nopu")[0].split("hEff_")[1] )])
        hlist[iplot].SetMarkerColor(markercolordir[(hlist[iplot].GetName().split("_pu")[0].split("hEff_")[1] if pued else hlist[iplot].GetName().split("_nopu")[0].split("hEff_")[1])])

        leg.AddEntry(hlist[iplot], namedir[(hlist[iplot].GetName().split("_pu")[0].split("hEff_")[1] if pued else hlist[iplot].GetName().split("_nopu")[0].split("hEff_")[1])], "P")
        #leg.AddEntry(hlist[iplot], namedir[(hlist[iplot].GetName().split("_pu")[0].split("hEff_")[1] if pued else hlist[iplot].GetName().split("_nopu")[0].split("hEff_")[1])] + " w/ failures" * ("nofail" not in hlist[iplot].GetName() and "noage" not in hlist[iplot].GetName() and "RPC" in hlist[iplot].GetName()), "P")

        #hlist[iplot].Draw("P,hist" + (iplot != 0) * "same")
        hlist[iplot].Draw("PE,hist" + (iplot != 0) * "same")

    leg.Draw("HISTP")

    textlist = []
    linelist = []
    for ich in range(4):
        textlist.append(r.TText(.17 + ich * 0.1975, 0.20, chambTag[ich]))
        textlist[-1].SetNDC(True)
        textlist[-1].Draw("same")
        if ich != 3:
            #if ich == 2:
            if ich == 2 or ich == 1:
                linelist.append(r.TLine(0.3075 + ich * 0.1975, 0.1, 0.3075 + ich * 0.1975, legylow))
                linelist[-1].SetNDC(True)
                linelist[-1].Draw("same")
                linelist.append(r.TLine(0.3075 + ich * 0.1975, legyhigh, 0.3075 + ich * 0.1975, 0.9))
                linelist[-1].SetNDC(True)
                linelist[-1].Draw("same")
            else:
                linelist.append(r.TLine(0.3075 + ich * 0.1975, 0.1, 0.3075 + ich * 0.1975, 0.9))
                linelist[-1].SetNDC(True)
                linelist[-1].Draw("same")

    #cmslat = r.TLatex()
    #cmslat.SetTextSize(0.03);
    #cmslat.DrawLatexNDC(0.11, 0.91, "#scale[1.5]{CMS}");
    #cmslat.Draw("same");

    #CMS_lumi.lumi_13TeV = ""
    #CMS_lumi.extraText  = 'Phase-2 Simulation'
    #CMS_lumi.cmsTextSize= 0.5
    #CMS_lumi.lumi_sqrtS = ''
    #CMS_lumi.CMS_lumi(r.gPad, 0, 0, 0.07)

    firsttex = r.TLatex()
    firsttex.SetTextSize(0.03)
    #firsttex.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Phase-2 Simulation")  ## La de Brieuc
    firsttex.DrawLatexNDC(0.11, 0.92, "#scale[1.5]{CMS} Phase-2 Simulation")
    firsttex.Draw("same");

    prelim_text = r.TLatex()
    prelim_text.SetTextSize(0.03)
    prelim_text.DrawLatexNDC(0.208, 0.865, "Preliminary")

    secondtext = r.TLatex()
    if pued: toDisplay = r.TString("14 TeV, 3000 fb^{-1}, 200 PU, p_{T}>20 GeV")
    else:    toDisplay = r.TString("14 TeV, 3000 fb^{-1}, 0 PU, p_{T}>20 GeV")

    #if pued: toDisplay = r.TString("14 TeV, 3000/1000 fb^{-1}, 200 PU, p_{T}>20 GeV")
    #else:    toDisplay = r.TString("14 TeV, 3000/1000 fb^{-1}, 0 PU, p_{T}>20 GeV")

    #secondtext.SetTextSize(0.035) ## La de Brieuc
    secondtext.SetTextSize(0.0275) # Good one
    #secondtext.SetTextSize(0.025)
    secondtext.SetTextAlign(31)
    #secondtext.DrawLatexNDC(0.90, 0.91, toDisplay.Data())  ## La de Brieuc
    secondtext.DrawLatexNDC(0.90, 0.92, toDisplay.Data())
    secondtext.Draw("same")

    #c.SetLogy()

    outputscaff = path + "/" + savescaffold.format(pu = (not pued) * "no" + "pu", qu = ("_" + qual)*(qual != ""), id =("_" + ind)*(ind != "")) + "_zoomout" * (zoom == "zoomout") + "_zoomin" * (zoom == "zoomin") + "_extremezoomin" * (zoom == "extremezoomin") + suffix + extrasuffix

    c.SaveAs(outputscaff + ".png")
    c.SaveAs(outputscaff + ".pdf")
    c.SaveAs(outputscaff + ".root")
    c.SaveAs(outputscaff + ".C")
    c.Close(); del c
    return


def producetheTDRplot(qual = "", pu = False, ind = "", zoom = True):
    print "\nBeginning plotting for pu", pu, "\n"
    listofplots = []

    #makeresplot(listofplots, "noage", "AM",     qual, "nopu" if not pu else "pu200", ind)
    #makeresplot(listofplots, "noage", "HB",     qual, "nopu" if not pu else "pu200", ind)
    #makeresplot(listofplots, "noage", "AM+RPC", qual, "nopu" if not pu else "pu200", ind)

    #makeresplot(listofplots, scenariolist[1],   "AM",     qual, "nopu" if not pu else "pu200", ind)
    #makeresplot(listofplots, scenariolist[1],   "HB",     qual, "nopu" if not pu else "pu200", ind)
    #makeresplot(listofplots, scenariolist[1],   "AM+RPC", qual, "nopu" if not pu else "pu200", ind)


    #makeresplot(listofplots, scenariolist[2],   "AM",     qual, "nopu" if not pu else "pu200", ind)
    #makeresplot(listofplots, scenariolist[2],   "HB",     qual, "nopu" if not pu else "pu200", ind)
    #makeresplot(listofplots, scenariolist[2],   "AM+RPC", qual, "nopu" if not pu else "pu200", ind)

    #makeresplot(listofplots, "noage", "AM",     qual, "nopu" if not pu else "pu200", ind)

    for scen in scenariolist:
        makeresplot(listofplots, scen, "AM",     qual, "nopu" if not pu else "pu200", ind)
        makeresplot(listofplots, scen, "HB",     qual, "nopu" if not pu else "pu200", ind)
        makeresplot(listofplots, scen, "AM+RPC", qual, "nopu" if not pu else "pu200", ind)

    print "\nCombining and saving\n"
    combineresplots(listofplots, qual, pu, ind, zoom)
    return


#def producetheSilviaplots():   ### NOT UPDATED
    #producetheTDRplot("", False, "ind0")
    #producetheTDRplot("", True,  "ind0")
    #producetheTDRplot("", False, "indmax2")
    #producetheTDRplot("", True,  "indmax2")
    #producetheTDRplot("higherthanfour", False)
    #producetheTDRplot("higherthanfour", True)
    #producetheTDRplot("higherthanfourvetoing", False)
    #producetheTDRplot("higherthanfourvetoing", True)
    #return


#### ACTUAL EXECUTION AND SMALL NICE CODE
if not os.path.isdir(path):
    print "Creating output folder..."
    os.system("mkdir " + path)
    os.system("cp " + path + "/../index.php " + path + "/")

#producetheTDRplot("")
#producetheTDRplot("", True)
#producetheTDRplot("nothreehits")
#producetheTDRplot("nothreehits", True)
producetheTDRplot("qualityOR")
producetheTDRplot("qualityOR", True)

#producetheTDRplot("",                  zoom = False)
#producetheTDRplot("", True,            zoom = False)
#producetheTDRplot("nothreehits",       zoom = False)
#producetheTDRplot("nothreehits", True, zoom = False)
producetheTDRplot("qualityOR",          zoom = "zoomout")
producetheTDRplot("qualityOR", True,    zoom = "zoomout")
producetheTDRplot("qualityOR",          zoom = "zoomin")
producetheTDRplot("qualityOR", True,    zoom = "zoomin")
producetheTDRplot("qualityOR",          zoom = "extremezoomin")
producetheTDRplot("qualityOR", True,    zoom = "extremezoomin")

#producetheSilviaplots()     ### NOT UPDATED
