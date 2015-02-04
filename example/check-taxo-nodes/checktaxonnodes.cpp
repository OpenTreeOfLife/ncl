#include <fstream>
#include "ncl/othelpers.h"
using namespace std;
bool gVerbose = false;
bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);
/* use some globals, because I'm being lazy... */
NxsSimpleTree * gRefTree = 0;
NxsSimpleTree * gTaxonTree = 0;
std::map<long, const NxsSimpleNode *> gOttID2RefNode;
std::map<const NxsSimpleNode *, std::string> gRefTipToName;
std::map<long, const NxsSimpleNode *> gOttID2TaxNode;
std::map<const NxsSimpleNode *, long> gTaxNode2ottID;
std::set<const NxsSimpleNode *> gSupportedNodes;
std::string gCurrentFilename;
std::string gCurrTmpFilepath;
//std::ostream * gCurrTmpOstream = 0L;


std::string getLeftmostDesName(const NxsSimpleNode *nd) {
	const std::string & name = nd->GetName();
	if (!name.empty()) {
		return name;
	}
	const unsigned outDegree = nd->GetChildren().size();
	if (outDegree == 0) {
		return gRefTipToName[nd];
	}
	return getLeftmostDesName(nd->GetChildren()[0]);
}

std::string getRightmostDesName(const NxsSimpleNode *nd) {
	const std::string & name = nd->GetName();
	if (!name.empty()) {
		return name;
	}
	const unsigned outDegree = nd->GetChildren().size();
	if (outDegree == 0) {
		return gRefTipToName[nd];
	}
	const unsigned lastInd = outDegree - 1;
	return getRightmostDesName(nd->GetChildren()[lastInd]);
}


void extendSupportedToRedundantNodes(const NxsSimpleTree * tree, std::set<const NxsSimpleNode *> & gSupportedNodes) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		std::vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree == 1 && gSupportedNodes.find(children[0]) != gSupportedNodes.end()) {
			gSupportedNodes.insert(nd);
		}
	}
}

bool singleDesSupportedOrNamed(const NxsSimpleNode *nd) {
	if (gSupportedNodes.find(nd) != gSupportedNodes.end()) {
		return true;
	}
	std::vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 1) {
		if (!nd->GetName().empty()) {
			return true;
		} else {
			return singleDesSupportedOrNamed(children[0]);
		}
	}
	return false;
}

std::map<const NxsSimpleNode *, std::set<long> > gRefNdp2mrca;
std::map<const NxsSimpleNode *, std::set<long> > gTaxNdp2mrca;
set<long> gRefLeafSet;
set<long> gTaxLeafSet;
map<const NxsSimpleNode *, long> gRefNamedNodes;


void processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, *nd);
		if (nd->IsTip()) {
			const unsigned ind = nd->GetTaxonIndex();
			assert(ind < tb->GetNumTaxonLabels());
			const std::string tn = tb->GetTaxonLabel(ind);
			gRefTipToName[nd] = tn;
			assert(ottID >= 0);
			gRefNdp2mrca[nd].insert(ottID);
			gRefLeafSet.insert(ottID);
		} else {
			useChildrenToFillMRCASet(nd, gRefNdp2mrca);
			if (ottID > 0) {
				gRefNamedNodes[nd] = ottID;
			}
		}
		if (ottID >= 0) {
			assert(gOttID2RefNode.find(ottID) == gOttID2RefNode.end());
			gOttID2RefNode[ottID] = nd;
		}
	}
}

void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(gOttID2TaxNode.find(ottID) == gOttID2TaxNode.end());
		if (nd->IsTip()) {
			assert(gOttID2RefNode.find(ottID) != gOttID2RefNode.end());
			gTaxLeafSet.insert(ottID);
			gTaxNdp2mrca[nd].insert(ottID);
		} else {
			useChildrenToFillMRCASet(nd, gTaxNdp2mrca);
		}
		gOttID2TaxNode[ottID] = *nIt;
		gTaxNode2ottID[*nIt] = ottID;
	}
	for (std::map<long, const NxsSimpleNode *>::const_iterator nit = gOttID2RefNode.begin(); nit != gOttID2RefNode.end(); ++nit) {
		assert(gOttID2TaxNode.find(nit->first) != gOttID2TaxNode.end());
	}
}


void writeSetDiff(std::ostream & out, const char *indent, const set<long> &fir, const char *firN, const set<long> & sec, const char *secN) {
		for (set<long>::const_iterator rIt = fir.begin(); rIt != fir.end(); ++rIt) {
			if (sec.find(*rIt) == sec.end()) {
				out << indent << "ott" << *rIt << " is in " << firN << " but not " << secN << "\n";
			}
		}
		for (set<long>::const_iterator rIt = sec.begin(); rIt != sec.end(); ++rIt) {
			if (fir.find(*rIt) == fir.end()) {
				out << indent << "ott" << *rIt << " is in " << secN << " but not " << firN << "\n";
			}
		}

}

bool isProperSubset(const set<long> & small, const set<long> & big) {
	if (big.size() <= small.size()) {
		return false;
	}
	for (set<long>::const_iterator rIt = small.begin(); rIt != small.end(); ++rIt) {
		if (big.find(*rIt) == big.end()) {
			return false;
		}
	}
	return true;
}

bool doCheckEquivalent(std::ostream &out, long ottID, const NxsSimpleNode * snode, std::map<const NxsSimpleNode *, std::set<long> > & srcLookup,
										  const NxsSimpleNode * tnode, std::map<const NxsSimpleNode *, std::set<long> > & taxLookup,
										  bool topLevel, bool climbSynth, bool climbTax) {
	std::map<const NxsSimpleNode *, std::set<long> >::const_iterator streeLSIt = srcLookup.find(snode);
	assert(streeLSIt != srcLookup.end());
	std::map<const NxsSimpleNode *, std::set<long> >::const_iterator taxtreeLSIt = taxLookup.find(tnode);
	assert(taxtreeLSIt != taxLookup.end());
	const std::set<long> & streeMRCA = streeLSIt->second;
	const std::set<long> & taxtreeMRCA = taxtreeLSIt->second;
	if (streeMRCA != taxtreeMRCA) {
		if (topLevel) {
			out << "ottID " << ottID << " incorrect:\n";
			writeSetDiff(out, "    ", streeMRCA, "synth", taxtreeMRCA, "taxonomy");
		}
		if (climbSynth && isProperSubset(streeMRCA, taxtreeMRCA)) {
			return doCheckEquivalent(out, ottID, snode->GetEdgeToParent().GetParent(), srcLookup, tnode, taxLookup, false, true, false);
		} else if (climbTax && isProperSubset(taxtreeMRCA, streeMRCA)) {
			return doCheckEquivalent(out, ottID, snode, srcLookup, tnode->GetEdgeToParent().GetParent(), taxLookup, false, false, true);
		} else {
			return false;
		}
	} else if (!topLevel) {
		out << "        Found identical leaf sets for the synthetic tree \"" << snode->GetName() << "\" and the taxonomic node \"" << tnode->GetName() << "\".\n";
	}
	return true;
}

void summarize(std::ostream & out) {
	for (map<const NxsSimpleNode *, long>::const_iterator rnit = gRefNamedNodes.begin(); rnit != gRefNamedNodes.end(); ++rnit) {
		const NxsSimpleNode * nd = rnit->first;
		const long ottID = rnit->second;
		std::map<long, const NxsSimpleNode *>::const_iterator tID2nd = gOttID2TaxNode.find(ottID);
		assert(tID2nd != gOttID2TaxNode.end());
		const NxsSimpleNode *taxNd = tID2nd->second;
		if (!doCheckEquivalent(out, ottID, nd, gRefNdp2mrca, taxNd, gTaxNdp2mrca, true, true, true)) {
			out << "        Could not find this set of leaves in the synth \"" << nd->GetName() <<"\" in any taxonomic node.\n";
		}
	}
	if (gTaxLeafSet != gRefLeafSet) {
		writeSetDiff(out, "", gRefLeafSet, "synth", gTaxLeafSet, "taxonomy");
	}
}

bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	static unsigned long gTreeCount = 0;
	const NxsTaxaBlockAPI * taxa = treesB->GetTaxaBlockPtr();
	gTreeCount++;
	if (gVerbose)
		std::cerr << "Read tree " <<  gTreeCount<< '\n';
	unsigned int nUnlabeledOutDegOne = 0;
	unsigned int nLabeledOutDegOne = 0;
	std::vector<std::string> parNames;
	NxsSimpleTree * nst = new NxsSimpleTree(ftd, 0.0, 0, true);
	//gExpanded.clear();
	//gTabooLeaf.clear();
	if (gRefTree == 0) {
		gRefTree = nst;
		processRefTree(taxa, nst);
	} else if (gTaxonTree == 0) {
		gTaxonTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		const char * msg = "Exepting only 2 files: the synthetic tree, and then the taxonomy\n";
		std::cerr << msg;
		throw NxsException(msg);
	}
	if (gRefTree != nst && gTaxonTree != nst) {
		delete nst;
	}
	return false;
}

void printHelp(ostream & out) {
	out << "chectaxonnodes tree.tre taxonomy.tre\n";
}

int main(int argc, char *argv[])
	{
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);

	bool readfile = false;
	bool el = false;
	bool depth = false;
	bool brief = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			printHelp(cout);
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'v')
			gVerbose = true;
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'f')
			{
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2)
				{
				std::string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after after -f\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		else
			{
			readfile = true;
			const std::string filepathstr(filepath);
			const size_t sp = filepathstr.find_last_of('/');
			if (sp == std::string::npos)
				{
				gCurrentFilename = filepathstr;
				}
			else
				{
				gCurrentFilename = filepathstr.substr(sp + 1);
				}
			gCurrTmpFilepath = std::string("tmp/") + gCurrentFilename;
			try {
				readFilepathAsNEXUS(filepath, f, newTreeHook, 0L);
			} catch (...) {
				//tostream.close();
				throw;
				}
			//gCurrTmpOstream = 0L;
			}
		}
	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
	summarize(std::cout);
	return 0;
	}

