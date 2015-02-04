#include <fstream>
#include "ncl/othelpers.h"
OTCLI gOTCli("findunsupportededgs: takes a tree, an optional file (.txt extension) of MRCA designators, a taxonomy, and the some # of newick files. It reports any edges in the first tree that have no support in the last trees. If the designators file is sent, then the nodes that are the MRCAs of the numbers on each line of that file will be treated as a priori problems to be checked. So you will get a report for those nodes whether or not they are unsupported.",
			 "findunsupportededgs synth.tre probs.txt taxonomy.tre one.tre two.tre three.tre");
using namespace std;
void processContent(PublicNexusReader & nexusReader, ostream *out);

bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);

/* use some globals, because I'm being lazy... */
NxsSimpleTree * gRefTree = 0;
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


void extendSupportedToRedundantNodes(const NxsSimpleTree * tree, set<const NxsSimpleNode *> & gSupportedNodes) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree == 1) {
			if (gSupportedNodes.find(children[0]) != gSupportedNodes.end() || children[0]->GetName().length() > 0) {
				if (gAPrioriProblemNode.find(nd) != gAPrioriProblemNode.end()) {
					assert(false); // shouldn't get out-degree one nodes w/ our designators
				}
				gSupportedNodes.insert(nd);
			}
		}
	}
}

bool singleDesSupportedOrNamed(const NxsSimpleNode *nd) {
	if (gSupportedNodes.find(nd) != gSupportedNodes.end()) {
		return true;
	}
	if (nd->GetOutDegree() == 1) {
		if (!nd->GetName().empty()) {
			return true;
		} else {
			return singleDesSupportedOrNamed(nd->GetFirstChild());
		}
	}
	return false;
}

int describeUnnamedUnsupported(ostream &out, const NxsSimpleTree * tree, const set<const NxsSimpleNode *> & supported) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin();
	int numUnsupported = 0;
	++nIt; //skip the root
	for (;nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree > 0 && supported.find(nd) == supported.end()) {
			if (IsRedundantNodeAroundTip(nd)) {
				//pass
			} else if (outDegree == 1 && gTrustNamedNodes && singleDesSupportedOrNamed(nd)) {
				//pass
			} else if (nd->GetName().length() == 0) { //assume that it is from the taxonomy
				if (gNoAprioriTests) {
					out << "Unsupported node ";
				} else {
					map<const NxsSimpleNode *, set<long> >::const_iterator gaIt = gAPrioriProblemNode.find(nd);
					if (gaIt == gAPrioriProblemNode.end()) {
						out << "Novel unsupported node ";
					} else {
						out << "Confirmation of unsupported node (designators =";
						writeOttSet(out, "", gaIt->second, " ");
						out << ") ";
					}
				}
				describeUnnamedNode(nd, out, 0, gRefTipToName, false);
				numUnsupported += 1;
			}
		}
	}
	return numUnsupported;
}


void summarize(ostream & out) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	extendSupportedToRedundantNodes(gRefTree, gSupportedNodes);
	int numUnsupported = describeUnnamedUnsupported(out, gRefTree, gSupportedNodes);
	for (map<const NxsSimpleNode *, set<long> >::const_iterator gaIt = gAPrioriProblemNode.begin(); gaIt != gAPrioriProblemNode.end(); ++gaIt) {
		if (gSupportedNodes.find(gaIt->first) != gSupportedNodes.end()) {
			out << "Claim of unsupported apparently refuted for designators: ";
			writeOttSet(out, "", gaIt->second, " ");
			out << ". See standard error stream for details.\n";
		}
	} 
	int supNonNamed = 0;
	int numSupportedInternals = 0;
	for (set<const NxsSimpleNode *>::const_iterator rIt = gSupportedNodes.begin(); rIt != gSupportedNodes.end(); ++rIt) {
		if (!(*rIt)->IsTip()) {
			numSupportedInternals += 1;
			if ((*rIt)->GetName().length() == 0) {
				supNonNamed += 1;
			}
		}
	}
	out << "\n\nFinal summary:\n";
	out << gRefTreeNumNamedInternalsNodes << " internal nodes were named in the reference tree. These were not rigorously checked against the taxonomy. They may not be detected as errors.\n";
	out << numSupportedInternals << " internal nodes where flagged as being supported by an input (including taxonomy).\n";
	int supNamed = numSupportedInternals - supNonNamed;
	out << "    " << supNamed << " of these were named (some of the support could just be the taxonomic expansion of tips).\n";
	out << "    " << supNonNamed << " of these were unnamed.\n";
	out << numUnsupported << " unsupported nodes.\n";
	out << endl;
	if (gExitCode < 0) {
		gExitCode -= numUnsupported;
	} else {
		gExitCode = numUnsupported;
	}
}




void processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, *nd);
		if (nd->IsTip()) {
			const unsigned ind = nd->GetTaxonIndex();
			assert(ind < tb->GetNumTaxonLabels());
			const string tn = tb->GetTaxonLabel(ind);
			gRefTipToName[nd] = tn;
		} else if (nd->GetName().length() > 0) {
			gRefTreeNumNamedInternalsNodes += 1;
		}
		if (ottID >= 0) {
			assert(gOttID2RefNode.find(ottID) == gOttID2RefNode.end());
			gOttID2RefNode[ottID] = nd;
		}
	}
}

void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		std::cerr << "nd->GetName()  = " << nd->GetName() << '\n';
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(gOttID2TaxNode.find(ottID) == gOttID2TaxNode.end());
		if (nd->IsTip()) {
			assert(gOttID2RefNode.find(ottID) != gOttID2RefNode.end());
			gTaxLeafOTTIDs.insert(ottID);
		}
		gOttID2TaxNode[ottID] = *nIt;
		gTaxNode2ottID[*nIt] = ottID;
	}
	for (map<long, const NxsSimpleNode *>::const_iterator nit = gOttID2RefNode.begin(); nit != gOttID2RefNode.end(); ++nit) {
		assert(gOttID2TaxNode.find(nit->first) != gOttID2TaxNode.end());
	}
}



const NxsSimpleNode * findMRCA(const map<long, const NxsSimpleNode *> & ref,
	                           const map<long, const NxsSimpleNode *> &taxonomy, long ottID) {
	set<long> tipOTTIDs;
	fillTipOTTIDs(taxonomy, ottID, tipOTTIDs, gTaxNode2ottID);
	const unsigned nTips = tipOTTIDs.size();
	if (nTips < 2) {
		cerr << "findMRCA called on " << ottID << '\n';
		assert(false);
	}
	gNonMono[ottID] = tipOTTIDs;
	return findMRCAFromIDSet(ref, tipOTTIDs, ottID);
}


void markPathToRoot(map<const NxsSimpleNode *, set<long> > &n2m, long ottID) {
	const NxsSimpleNode * nd = 0;
	map<long, const NxsSimpleNode *>::const_iterator fIt = gOttID2RefNode.find(ottID);
	if (fIt != gOttID2RefNode.end()) {
		nd = fIt->second;
	} else {
		assert(false);
		nd = findMRCA(gOttID2RefNode, gOttID2TaxNode, ottID);
	}
	assert(nd != 0);
	while (nd != 0) {
		n2m[nd].insert(ottID);
		nd = nd->GetEdgeToParentRef().GetParent();
	}
}


void recordSupportedNodes(const map<const NxsSimpleNode *, set<long> > & refNdp2mrca,
						  const map<set<long>, const NxsSimpleNode *> & sourceClades,
						  set<const NxsSimpleNode *> & supportedNodes,
						  const set<long> & leafSet,
						  const map<const NxsSimpleNode *, set<long> > & srcNdp2mrca) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	if (false) {//debugging
		cerr << "sourceClades:\n";
		for (map<set<long>, const NxsSimpleNode *>::const_iterator scit = sourceClades.begin(); scit != sourceClades.end(); ++scit) {
			cerr << "  clade: ";
			writeOttSet(cerr, "", scit->first, " ");
			cerr << "\n";
		}
	}

	for (map<const NxsSimpleNode *, set<long> >::const_iterator nsIt = refNdp2mrca.begin(); nsIt != refNdp2mrca.end(); ++nsIt) {
		const NxsSimpleNode * nd = nsIt->first;
		const bool printDB = false; //gAPrioriProblemNode.find(nd) != gAPrioriProblemNode.end();
		const NxsSimpleNode * par = nd->GetEdgeToParentRef().GetParent();
		if (par != 0) {
			const set<long> & nm = nsIt->second;
			const NxsSimpleNode * firstBranchingAnc = findFirstBranchingAnc(nd);
			map<const NxsSimpleNode *, set<long> >::const_iterator ancIt = refNdp2mrca.find(firstBranchingAnc);
			assert(ancIt != refNdp2mrca.end());
			const set<long> & anm = ancIt->second;
			if (printDB) { //debugging
				cerr << "DEBUGGING refTree node " << (long)nd << ": ";
				writeOttSet(cerr, "", nm, " ");
				cerr << "\n";
				cerr << "refTree par " << (long)firstBranchingAnc << ": ";
				writeOttSet(cerr, "", anm, " ");
				cerr << "\n";
			}
			if (isTheMrcaInThisLeafSet(nd, refNdp2mrca)) {
				if (anm != nm) {
					map<set<long>, const NxsSimpleNode *>::const_iterator scIt = sourceClades.find(nm);
					if (scIt != sourceClades.end()) {
						if (gAPrioriProblemNode.find(nd) != gAPrioriProblemNode.end()) {
							ostream & out = cerr;
							map<const NxsSimpleNode *, set<long> >::const_iterator apIt = gAPrioriProblemNode.find(nd);
							out << "ERROR!: a priori unsupported node found. Designators were ";
							writeOttSet(out, "", apIt->second, " ");
							out << ". A node was found, which (when pruned to the leaf set of an input tree) contained:\n";
							writeOttSet(out, "    ", nm, " ");
							out << "\nThe subtree from the source was: ";
							const NxsSimpleNode * srcNd = scIt->second;
							writeSubtreeNewickOTTIDs(out, srcNd, srcNdp2mrca);
							gExitCode = -1;
						}
						supportedNodes.insert(nd);
						if (printDB) { //debugging
							cerr <<"SUPPORTED!\n";
						}
					} else {
						if (printDB) { //debugging
							cerr <<"UNsupported. node mrca set not in source tree\n";
						}
					} 
				} else {
					if (printDB) { //debugging
						cerr <<"UNsupported. node mrca set == first branching ancestors mrca set on this leaf set\n";
					}
				}
			} else {
				if (printDB) { //debugging
					cerr <<"UNsupported. not a mrcaInThisLeafSet\n";
				}
			}
		}
	}
}


void expandOTTInternalsWhichAreLeaves(const NxsTaxaBlockAPI * tb, NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	map<NxsSimpleNode *, set<long> > replaceNodes;
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode & nd = **nIt;
		const unsigned outDegree = nd.GetOutDegree();
		if (outDegree == 0) {
			long ottID = getOTTIndex(tb, **nIt);
			assert(ottID >= 0);
			assert(gOttID2TaxNode.find(ottID) != gOttID2TaxNode.end());
			if (gTaxLeafOTTIDs.find(ottID) == gTaxLeafOTTIDs.end()) {
				set<long> leafSet;
				fillTipOTTIDs(gOttID2TaxNode, ottID, leafSet, gTaxNode2ottID);
				replaceNodes[const_cast<NxsSimpleNode *>(&nd)] = leafSet;
			}
		}
	}
	for (map<NxsSimpleNode *, set<long> >::const_iterator rIt = replaceNodes.begin(); rIt != replaceNodes.end(); ++rIt) {
		NxsSimpleNode * oldNode = rIt->first;
		const set<long> & leafSet = rIt->second;
		assert(leafSet.size() > 0);
		oldNode->SetTaxonIndex(UINT_MAX); // make this no longer appear to be a tip
		for (set<long>::const_iterator lsIt = leafSet.begin(); lsIt != leafSet.end(); ++lsIt) {
			NxsSimpleNode *newNode =  tree->AllocNewNode(oldNode);
			oldNode->AddChild(newNode);
			gExpanded[newNode] = *lsIt;
			assert(gTabooLeaf.find(*lsIt) == gTabooLeaf.end());
			gTabooLeaf[*lsIt] = newNode;
		}
	}
}

void processSourceTree(const NxsTaxaBlockAPI * tb, NxsSimpleTree * tree) {
	expandOTTInternalsWhichAreLeaves(tb, tree);
	map<const NxsSimpleNode *, set<long> > ndp2mrca;
	map<const NxsSimpleNode *, set<long> > refNdp2mrca;
	set<long> leafSet;

	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode & nd = **nIt;
		vector<NxsSimpleNode *> children = nd.GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree > 0) {
			set<long> & mrca = ndp2mrca[&nd];
			for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
				const NxsSimpleNode * c = *cIt;
				assert(ndp2mrca.find(c) != ndp2mrca.end());
				set<long> & csl = ndp2mrca[c];
				mrca.insert(csl.begin(), csl.end());
			}
		} else if (outDegree == 0) {
			long ottID = getOTTIndex(tb, **nIt);
			map<long, const NxsSimpleNode *>::const_iterator tlIt = gTabooLeaf.find(ottID);
			if (tlIt != gTabooLeaf.end()) {
				assert(tlIt->second == *nIt);
			}
			assert(ottID >= 0);
			assert(gOttID2TaxNode.find(ottID) != gOttID2TaxNode.end());
			ndp2mrca[&nd].insert(ottID);
			markPathToRoot(refNdp2mrca, ottID);
			leafSet.insert(ottID);
		}
	}
	map<set<long>, const NxsSimpleNode *> sourceClades;
	for (map<const NxsSimpleNode *, set<long> >::const_iterator nsIt = ndp2mrca.begin(); nsIt != ndp2mrca.end(); ++nsIt) {
		const NxsSimpleNode * nd = nsIt->first;
		const NxsSimpleNode * par = nd->GetEdgeToParentRef().GetParent();
		if (par != 0 && !nd->IsTip()) {
			const set<long> & nm = nsIt->second;
			sourceClades[nm] = nsIt->first;
		}
	}
	recordSupportedNodes(refNdp2mrca, sourceClades, gSupportedNodes, leafSet, ndp2mrca);
	if (gCurrTmpOstream != 0) {
		*gCurrTmpOstream << "#NEXUS\nBEGIN TREES;\n";

		*gCurrTmpOstream << "   Tree pruned_synth = [&R] ";
		writeNewickOTTIDs(*gCurrTmpOstream, gRefTree, refNdp2mrca);
		
		*gCurrTmpOstream << "   Tree input = [&R] ";
		writeNewickOTTIDs(*gCurrTmpOstream, tree, ndp2mrca);
		
		*gCurrTmpOstream << "END;\n";
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
	if (gRefTree == 0) {
		gRefTree = nst;
		processRefTree(taxa, nst);
	} else if (gTaxonTree == 0) {
		gTaxonTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		processSourceTree(taxa, nst);
	}
	if (gRefTree != nst && gTaxonTree != nst) {
		delete nst;
	}
	return false;
}

void markSuspectNode(const set<long> & designators) {
	const NxsSimpleNode * mrca = findMRCAFromIDSet(gOttID2RefNode, designators, -1);
	assert(mrca->GetName().length() == 0);
	gAPrioriProblemNode[mrca] = designators;
}


void parseAndProcessMRCADesignatorsFile(string filepath) {
	if (gRefTree == 0L || gTaxonTree != 0L) {
		cerr << "gRefTree" << (long )gRefTree << "\n";
		cerr << "gTaxonTree" << (long )gTaxonTree << "\n";
		cerr << "\nDesignators file must come after the full tree estimate, but before the taxonomy in the argument length\n";
		throw exception();
	}
	ifstream inpf(filepath.c_str());
	assert(inpf.good());
	string line;
	while (getline(inpf, line)) {
		string stripped = NxsString::strip_surrounding_whitespace(line);
		if (!stripped.empty()) {
			list<string> words;
			NxsString::split(stripped, &words);
			if (words.size() < 2) {
				cerr << "Expecting >1 designator. Found: " << line << "\n";
				throw exception();
			}
			set<long> designators;
			for (list<string>::const_iterator dIt = words.begin(); dIt != words.end(); ++dIt) {
				long d;
				if (!NxsString::to_long(dIt->c_str(), &d)) {
					cerr << "Expecting numeric designator. Found: " << line << "\n";
					throw exception();
				}
				designators.insert(d);
			}
			markSuspectNode(designators);
		}
	}
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	inpf.close();
}

int main(int argc, char *argv[])
	{
	std::vector<std::string> args;
	if (!gOTCli.parseArgs(argc, argv, args)) {
		return 1;
	}
	const char * tooFewArgsMsg = "Expecting a complete tree file, a taxonomy tree, and at least one input tree file.\n";
	if (args.size() < 3) {
		cerr << tooFewArgsMsg;
		return 2;
	}
	try{
		gOTCli.readFilepath(args[0], newTreeHook);
		int n = 1;
		if (gOTCli.isDotTxtFile(args[n])) {
			gNoAprioriTests = false;
			parseAndProcessMRCADesignatorsFile(args[n]);
			n++;
		}
		if (n >= args.size() - 1) { // need to read taxonomy and then at least 1 tree
			cerr << tooFewArgsMsg;
			return 2;
		}
		for(; n < args.size(); ++n) {
			gOTCli.readFilepath(args[n], newTreeHook);
		}
		
	} catch (NxsException &x) {
		std::cerr << x.what() << "\n";
		return 1;
	}
	if (gOTCli.exitCode == 0) {
		summarize(std::cout);
	}
	return gOTCli.exitCode;
}
