# -*- coding: utf-8 -*-
import ROOT as r
from copy import deepcopy
import CMS_lumi, os
r.gROOT.SetBatch(True)

#path = "~sscruz/www/DT_TDR/2019_12_09_plots_eff_withHBaged_noquality/"
#path = "~vrbouza/www/Miscelánea/2019_09_23_plots_eff_shiftsoff_paSilvia"
path = "~vrbouza/www/Miscelánea/2019_10_04_plots_eff_shiftsoff"

plotscaffold = "hEff_{st}_{al}_{ty}"
savescaffold = "hEff_{pu}{qu}{id}"
chambTag = ["MB1", "MB2", "MB3", "MB4"]
suffix = ""

def makeresplot(hlist, aged, algo, qual = "", pued = False, ind = ""):
    print "\nObtaining intermediate plot for algo", algo, "which is", aged, "aged and considering", pued, "pile-up"
    print "The file that it's going to be opened is", "results_eff_" + ((not pued) * "no") + "pu_" + (not aged) * "no" + "age_" + ("with" * ("RPC" in algo) + "no" * ("RPC" not in algo)) + "rpc_" + qual + ind + suffix + ".root"

    res = r.TFile.Open("results_eff_" + ((not pued) * "no") + "pu_" + (not aged) * "no" + "age_" + ("with" * ("RPC" in algo) + "no" * ("RPC" not in algo)) + "rpc_{q}".format(q = qual) + "_{id}".format(id = ind) * (ind != "") + suffix + ".root")

    print "Then, we're gonna get the histograms labeled as follows:", plotscaffold.format(al = algo.replace("+RPC", ""), st = chambTag[0], ty = "matched")

    hmatched = [res.Get(plotscaffold.format(al = algo.replace("+RPC", ""), st = chambTag[ich], ty = "matched")) for ich in range(4)]
    htotal   = [res.Get(plotscaffold.format(al = algo.replace("+RPC", ""), st = chambTag[ich], ty = "total")) for ich in range(4)]

    print "And finally the histogram that is going to be saved will have this name:", "hEff_{al}_{ag}_{pu}".format(al = algo, ag = (not aged) * "no" + "age", pu = ((not pued) * "no") + "pu") + "_{q}".format(q = qual) * (qual != "") + "_{id}".format(id = ind)* (ind != "")

    resplot = r.TH1D("hEff_{al}_{ag}_{pu}".format(al = algo, ag = (not aged) * "no" + "age", pu = ((not pued) * "no") + "pu") + "_{q}".format(q = qual) * (qual != "") + "_{id}".format(id = ind)* (ind != ""), "", 20, -0.5, 19.5)
    
    ibin = 1
    for ich in range(4):
        for iwh in range(1, 6):
            #print "st", ich, "wh", iwh, "valmatched:", hmatched[ich].GetBinContent(iwh), "valtotal:", htotal[ich].GetBinContent(iwh)
            effval = hmatched[ich].GetBinContent(iwh) / htotal[ich].GetBinContent(iwh)

            if effval > 1: print "WARNING: eff. for chamber", ich, "and wheel", iwh, "is", effval, ". Fixing to 1."

            resplot.SetBinContent(ibin, effval * (effval < 1) + 1 * (effval >= 1))
            eff = r.TEfficiency('kk','',1,-0.5,0.5)
            eff.SetTotalEvents(1, int(htotal[ich].GetBinContent(iwh)))
            eff.SetPassedEvents(1,int(hmatched[ich].GetBinContent(iwh)))
            if (eff.GetEfficiencyErrorLow(1)-eff.GetEfficiencyErrorUp(1)) > 0.05: print 'warning, bin asymmetric'
            resplot.SetBinError( ibin, max(eff.GetEfficiencyErrorLow(1),eff.GetEfficiencyErrorUp(1)))
            del eff
            ibin += 1

    hlist.append(deepcopy(resplot))
    res.Close(); del hmatched, htotal, res, resplot
    return


#lowlimityaxis  = 0.2
lowlimityaxis  = 0.8
highlimityaxis = 1
markersize     = 1
yaxistitle     = "Efficiency (adim.)"
yaxistitleoffset= 1.5
xaxistitle     = "Wheel"
legxlow        = 0.3075 + 2 * 0.1975
legylow        = 0.3
legxhigh       = 0.9
legyhigh       = 0.5

markertypedir  = {}
markertypedir["AM_age"]       = 24
markertypedir["AM_noage"]     = 20
#markertypedir["AM+RPC_age"]   = 29
#markertypedir["AM+RPC_noage"] = 29
markertypedir["AM+RPC_age"]   = 28
markertypedir["AM+RPC_noage"] = 34
markertypedir["HB_noage"]     = 22
markertypedir["HB_age"]       = 22

markercolordir  = {}
markercolordir["AM_age"]       = 2
markercolordir["AM+RPC_age"]   = 2
markercolordir["AM+RPC_noage"] = 4
markercolordir["HB_noage"]     = 4
markercolordir["AM_noage"]     = 4
markercolordir["HB_age"]       = 2

namedir  = {}
namedir["AM_age"]       = "AM w/ ageing"
namedir["AM+RPC_age"]   = "AM w/ RPC w/ ageing"
namedir["AM+RPC_noage"] = "AM w/ RPC"
namedir["HB_noage"]     = "HB"
namedir["AM_noage"]     = "AM"
namedir["HB_age"]       = "HB w/ ageing"

def combineresplots(hlist, qual = "", pued = False, ind = ""):
    print "Combining list of plots that has", pued, "pile-up,", qual, "quality and", ind, "index."
    if len(hlist) == 0: raise RuntimeError("Empty list of plots")

    c   = r.TCanvas("c", "c", 800, 800)
    c.SetLeftMargin(0.11)
    c.SetGrid()

    leg = r.TLegend(legxlow, legylow, legxhigh, legyhigh)
    hlist[0].SetStats(False)
    #hlist[0].SetTitle("L1 DT Phase 2 algorithm efficiency comparison")
    hlist[0].GetYaxis().SetRangeUser(lowlimityaxis, highlimityaxis)
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
        hlist[iplot].SetMarkerStyle(markertypedir[( hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2] )])
        hlist[iplot].SetMarkerColor(markercolordir[hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2]])
        #leg.AddEntry(hlist[iplot], namedir[hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2]], "P")
        leg.AddEntry(hlist[iplot], namedir[hlist[iplot].GetName().split("_")[1] + "_" + hlist[iplot].GetName().split("_")[2]])
        hlist[iplot].Draw("P,hist" + (iplot != 0) * "same")

    leg.Draw("HISTP")

    textlist = []
    linelist = []
    for ich in range(4):
        textlist.append(r.TText(.17 + ich * 0.1975, 0.20, chambTag[ich]))
        textlist[-1].SetNDC(True)
        textlist[-1].Draw("same")
        if ich != 3:
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
    firsttex.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Phase-2 Simulation")
    firsttex.Draw("same");

    secondtext = r.TLatex()
    toDisplay  = r.TString("14 TeV, 3000 fb^{-1}, 200 PU")
    secondtext.SetTextSize(0.035)
    secondtext.SetTextAlign(31)
    secondtext.DrawLatexNDC(0.90, 0.91, toDisplay.Data())
    secondtext.Draw("same")

    #c.SetLogy()

    outputscaff = path + "/" + savescaffold.format(pu = (not pued) * "no" + "pu", qu = ("_" + qual)*(qual != ""), id =("_" + ind)*(ind != ""))

    c.SaveAs(outputscaff + ".png")
    c.SaveAs(outputscaff + ".pdf")
    c.SaveAs(outputscaff + ".root")
    c.SaveAs(outputscaff + ".C")
    c.Close(); del c
    return


def producetheTDRplot(qual = "", pu = False, ind = ""):
    print "\nBeginning plotting for pu", pu, "\n"
    listofplots = []
    makeresplot(listofplots, False, "AM",     qual, pu, ind)
    makeresplot(listofplots, True,  "AM",     qual, pu, ind)
    makeresplot(listofplots, False, "HB",     qual, pu, ind)
    makeresplot(listofplots, True,  "HB",     qual, pu, ind)
    makeresplot(listofplots, False, "AM+RPC", qual, pu, ind)
    makeresplot(listofplots, True,  "AM+RPC", qual, pu, ind)

    print "\nCombining and saving\n"
    combineresplots(listofplots, qual, pu, ind)


def producetheSilviaplots():
    producetheTDRplot("", False, "ind0")
    producetheTDRplot("", True,  "ind0")
    producetheTDRplot("", False, "indmax2")
    producetheTDRplot("", True,  "indmax2")
    producetheTDRplot("higherthanfour", False)
    producetheTDRplot("higherthanfour", True)
    producetheTDRplot("higherthanfourvetoing", False)
    producetheTDRplot("higherthanfourvetoing", True)


#### ACTUAL EXECUTION AND SMALL NICE CODE

if not os.path.isdir(path):
    print "Creating output folder..."
    os.system("mkdir " + path)
    os.system("cp " + path + "/../index.php " + path + "/")

producetheTDRplot("")
producetheTDRplot("", True)
producetheTDRplot("nothreehits")
producetheTDRplot("nothreehits", True)
producetheTDRplot("qualityOR")
producetheTDRplot("qualityOR", True)

#producetheSilviaplots()
