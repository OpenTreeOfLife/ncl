#include "ncl/othelpers.h"
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
using namespace std;
long gStrictLevel = 2;
////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
////////////////////////////////////////////////////////////////////////////////
void processFilepath(
	const char * filename, // name of the file to be read
	ostream *out, // output stream to use (NULL for no output). Not that cerr is used to report errors.
	MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
	ProcessedTreeValidationFunction func, /*!< your pointer to your callback function */
	void * blob) {
	assert(filename);
	try {
		MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
		NxsCharactersBlock * charsB = nexusReader.GetCharactersBlockTemplate();
		NxsDataBlock * dataB = nexusReader.GetDataBlockTemplate();
		charsB->SetAllowAugmentingOfSequenceSymbols(true);
		dataB->SetAllowAugmentingOfSequenceSymbols(true);
		NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
		assert(treesB);
		if (gStrictLevel < 2)
			treesB->SetAllowImplicitNames(true);
		treesB->setValidationCallbacks(func, blob);
		if (gStrictLevel < 2) {
			NxsStoreTokensBlockReader *storerB =  nexusReader.GetUnknownBlockTemplate();
			assert(storerB);
			storerB->SetTolerateEOFInBlock(true);
		}
		cerr << "Executing" <<endl;
		try {
			nexusReader.ReadFilepath(filename, fmt);
		} catch(...) {
			nexusReader.DeleteBlocksFromFactories();
			throw;
		}
		nexusReader.DeleteBlocksFromFactories();
	} catch (const NxsException &x) {
		cerr << "Error:\n " << x.msg << endl;
		if (x.line >=0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
	}
}

int readFilepathAsNEXUS(const char *filename,
						MultiFormatReader::DataFormatType fmt,
						ProcessedTreeValidationFunction func, 
						void * blob) {
	cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = 0L;
		outStream = &cout;
		processFilepath(filename, outStream, fmt, func, blob);
	} catch (...) {
		cerr << "Parsing of " << filename << " failed (with an exception)" << endl;
		return 1;
	}
	return 0;
}


/*! \returns 0 on success*/
int readFilesListedInFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt,
						ProcessedTreeValidationFunction func, /*!< your pointer to your callback function */
	  				  	void * blob)
	{
	ifstream masterStream(masterFilepath, ios::binary);
	if (masterStream.bad())
		{
		cerr << "Could not open " << masterFilepath << "." << endl;
		exit(3);
		}
	char filename[1024];
	while ((!masterStream.eof())  && masterStream.good())
		{
		masterStream.getline(filename, 1024);
		if (strlen(filename) > 0 && filename[0] != '#')
			{
			int rc = readFilepathAsNEXUS(filename, fmt, func, blob);
			if (rc != 0)
				return rc;
			}
		}
	return 0;
	}


void OTCLI::printHelp(std::ostream & out) {
	out << this->titleStr << "\n";
	out << "\nThe most common usage is simply:\n";
	out << this->usageStr << "\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -v verbose output\n\n";
	out << "    -f<format> specifies the input file format expected:\n";
	out << "            -fnexus     NEXUS (this is also the default)\n";
	out << "            -frelaxedphyliptree  newick(this is also the default)\n";
	out << "        The complete list of format names that can follow the -f flag is:\n";
	vector<string> fmtNames =  MultiFormatReader::getFormatNames();
	for (vector<string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n) {
		out << "            "<< *n << "\n";
	}
}

bool OTCLI::parseArgs(int argc, char *argv[], std::vector<std::string> & args) {
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h') {
			this->printHelp(cout);
			this->exitCode = 1;
			return false;
		} else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'v') {
			this->verbose = true;
		} else if (slen > 1 && filepath[0] == '-' && filepath[1] == 's') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &(this->strictLevel)))) {
				cerr << "Expecting an integer after -s\n" << endl;
				this->printHelp(cerr);
				this->exitCode = 2;
				return false;
			}
		} else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'f') {
			fmt = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2) {
				string fmtName(filepath + 2, slen - 2);
				fmt =  MultiFormatReader::formatNameToCode(fmtName);
			}
			if (fmt == MultiFormatReader::UNSUPPORTED_FORMAT) {
				cerr << "Expecting a format after after -f\n" << endl;
				printHelp(cerr);
				this->exitCode = 2;
				return false;
			}
		} else {
			const string filepathstr(filepath);
			args.push_back(filepathstr);
		}
	}
	this->exitCode = 0;
	return true;
}

bool OTCLI::isDotTxtFile(const std::string &fp) {
	const size_t fnl = fp.length();
	return (fp.substr(fnl - 4) == string(".txt"));
}

int OTCLI::readFilepath(const std::string &fp, ProcessedTreeValidationFunction func, void * blob) {
	const size_t sp = fp.find_last_of('/');
	if (sp == std::string::npos) {
		this->currentFilename = fp;
	} else {
		this->currentFilename = fp.substr(sp + 1);
	}
	this->currTmpFilepath = std::string("tmp/") + this->currentFilename;
	this->currReadingDotTxtFile = this->isDotTxtFile(this->currentFilename);
	return readFilepathAsNEXUS(fp.c_str(), this->fmt, func, blob);
}

//@recursive!
void fillTipOTTIDs(const map<long, const NxsSimpleNode *> &taxonomy,
					long ottID,
					set<long> & tipOTTIDs,
					map<const NxsSimpleNode*, long> & taxNode2ottID) {
	map<long, const NxsSimpleNode *>::const_iterator tnIt = taxonomy.find(ottID);
	assert(tnIt != taxonomy.end());
	const NxsSimpleNode * tn = tnIt->second;
	vector<NxsSimpleNode *> children = tn->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 0) {
		tipOTTIDs.insert(taxNode2ottID[tn]);
	} else {
		for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			long ct = taxNode2ottID[*cIt];
			fillTipOTTIDs(taxonomy, ct, tipOTTIDs, taxNode2ottID);
		}
	}
}

bool isTheMrcaInThisLeafSet(const NxsSimpleNode * nd,
					   const map<const NxsSimpleNode *, set<long> > & refNdp2mrcaThisLeafSet) {
	vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree < 2) {
		return false;
	}
	bool foundFirstInf = false;
	for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
		const NxsSimpleNode * c = *cIt;
		const map<const NxsSimpleNode *, set<long> >::const_iterator rmIt = refNdp2mrcaThisLeafSet.find(c);
		if (rmIt != refNdp2mrcaThisLeafSet.end()) {
			if (foundFirstInf) {
				return true;
			}
			foundFirstInf = true;
		}
	}
	return false;
}
