#include <fstream>
#include "ncl/othelpers.h"
OTCLI gOTCli("tree: takes a taxonomy, and the some # of newick files. It reports whether trees have higher taxa as tips. If the designators file is sent, then the nodes that are the MRCAs of the numbers on each line of that file will be treated as a priori problems to be checked. So you will get a report for those nodes whether or not they are unsupported.",
			 "findhighertaxontips taxonomy.tre one.tre two.tre three.tre");
using namespace std;
void processContent(PublicNexusReader & nexusReader, ostream *out);

bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);

/* use some globals, because I'm being lazy... */
NxsSimpleTree * gTaxonTree = 0;
map<long, const NxsSimpleNode *> gOttID2RefNode;
map<const NxsSimpleNode *, string> gRefTipToName;
map<long, const NxsSimpleNode *> gOttID2TaxNode;
map<const NxsSimpleNode *, long> gTaxNode2ottID;
set<const NxsSimpleNode *> gSupportedNodes;
string gCurrentFilename;
string gCurrTmpFilepath;
ostream * gCurrTmpOstream = 0L;
bool gReadingTxtFile = false;
map<long, set<long> > gNonMono;
const bool gTrustNamedNodes = true;
map<const NxsSimpleNode *, long> gExpanded;
map<long, const NxsSimpleNode *> gTabooLeaf;
set<long> gTaxLeafOTTIDs;
map<const NxsSimpleNode *, set<long> > gAPrioriProblemNode;
bool gNoAprioriTests = true;
int gRefTreeNumNamedInternalsNodes = 0;
int gExitCode = 0;


void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		/*std::cerr << "nd->GetName()  = " << nd->GetName() << '\n';*/
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(gOttID2TaxNode.find(ottID) == gOttID2TaxNode.end());
		if (nd->IsTip()) {
			assert(ottID >= 0);
			gOttID2RefNode[ottID] = nd;
			gTaxLeafOTTIDs.insert(ottID);
		}
		gOttID2TaxNode[ottID] = *nIt;
		gTaxNode2ottID[*nIt] = ottID;
	}
	for (map<long, const NxsSimpleNode *>::const_iterator nit = gOttID2RefNode.begin(); nit != gOttID2RefNode.end(); ++nit) {
		assert(gOttID2TaxNode.find(nit->first) != gOttID2TaxNode.end());
	}
}



void treehasleaveswhichareOTTInternals(const NxsTaxaBlockAPI * tb, NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	map<NxsSimpleNode *, set<long> > replaceNodes;
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode & nd = **nIt;
		const unsigned outDegree = nd.GetOutDegree();
		if (outDegree == 0) {
			long ottID = getOTTIndex(tb, **nIt);
			assert (ottID >= 0);
   			assert(gOttID2TaxNode.find(ottID) != gOttID2TaxNode.end());
			if (gTaxLeafOTTIDs.find(ottID) == gTaxLeafOTTIDs.end()) {
            		std::cout << "Higher taxon " << ottID << " is tip\n";
		    }
			}
		}
	}



bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	static unsigned long gTreeCount = 0;
	const NxsTaxaBlockAPI * taxa = treesB->GetTaxaBlockPtr();
	gTreeCount++;
	if (gOTCli.verbose) {
		cerr << "Read tree " <<  gTreeCount<< '\n';
	}
	unsigned int nUnlabeledOutDegOne = 0;
	unsigned int nLabeledOutDegOne = 0;
	vector<string> parNames;
	NxsSimpleTree * nst = new NxsSimpleTree(ftd, 0.0, 0, true);
	gExpanded.clear();
	gTabooLeaf.clear();
	if (gTaxonTree == 0) {
		gTaxonTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		treehasleaveswhichareOTTInternals(taxa, nst);
	}
	return false;
}

int main(int argc, char *argv[])
	{
	std::vector<std::string> args;
	if (!gOTCli.parseArgs(argc, argv, args)) {
		return 1;
	}
	const char * tooFewArgsMsg = "Expecting a a taxonomy tree, and at least one input tree file.\n";
	if (args.size() < 2) {
		cerr << tooFewArgsMsg;
		return 2;
	}
	try{
		gOTCli.readFilepath(args[0], newTreeHook);
		int n = 1;
		if (n >= args.size()) { // need to read taxonomy and then at least 1 tree
			cerr << tooFewArgsMsg;
			return 2;
		}
		for(; n < args.size(); ++n) {
			gOTCli.readFilepath(args[n], newTreeHook);
			std::cout << args[n] << " \n";
		}

	} catch (NxsException &x) {
		std::cerr << x.what() << "\n";
		return 1;
	}
	return gOTCli.exitCode;
}
