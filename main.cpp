#include <set>
#include <iostream>
#include <time.h>
#include <map>
#include <string>
#include <fstream>


#include "node.h"
#include "newicklex.h"
#include "treeutil.h"
#include "ewah/ewah.h"

#include "BipartiteMWIS.h"

using namespace std;


typedef ewah::EWAHBoolArray<uint32_t> bitmap;





struct NodeInfo {
	int id;
	bitmap clade;
	bitmap clade_comp;	//complement of the clade
};


struct TreePairInfo {
	Node* t1;
	Node* t2;
	map<string, int> label_to_leafid;
	map<Node*, NodeInfo> infos;


	bitmap _all_ones;	//temp variable

	int cpt_id, leaf_id;
	int nb_leaves;

	TreePairInfo(Node* t1, Node* t2) : t1(t1), t2(t2) {
		cpt_id = 1;
		leaf_id = 1;

		
		//count leaves, losing time...
		nb_leaves = 0;
		for (auto it = t1->begin(); it != t1->end(); ++it) {
			if ((*it)->is_leaf()) {
				nb_leaves++;
				_all_ones.set(nb_leaves);
			}
		}
		

		preprocess_tree_rec(t1, true);
		
		preprocess_tree_rec(t2, false);
	}

	void preprocess_tree_rec(Node* v, bool is_tree1) {
		
		infos[v].id = cpt_id;
		cpt_id++;
		
		if (v->is_leaf()) {
			
			if (is_tree1) {
				infos[v].id = leaf_id;
				
				infos[v].clade.set(leaf_id);
				infos[v].clade_comp = _all_ones.logicalandnot(infos[v].clade);

				label_to_leafid[v->label] = leaf_id;
				leaf_id++;
			}
			else {
				infos[v].id = label_to_leafid[v->label];
				infos[v].clade.set(infos[v].id);
				infos[v].clade_comp = _all_ones.logicalandnot(infos[v].clade);
			}
			
		}
		else {
			for (int i = 0; i < v->get_nb_children(); ++i) {
				preprocess_tree_rec(v->get_child(i), is_tree1);

				if (i == 0)
					infos[v].clade = infos[v->get_child(i)].clade;
				else 
					infos[v].clade = infos[v].clade | infos[v->get_child(i)].clade;
			}
			infos[v].clade_comp = _all_ones.logicalandnot(infos[v].clade);
		}
	}



	map<Node*, vector<Node*>> get_imcompats() {

		map<Node*, vector<Node*>> ret;
		int cpt = 0;
		for (auto it1 = t1->begin(); it1 != t1->end(); ++it1) {

			cpt++;

			if (cpt % 100 == 0)
				cout << "cpt=" << cpt << endl;
			Node* v1 = (*it1);

			
			bitmap A = infos[v1].clade;
			bitmap B = infos[v1].clade_comp;

			for (auto it2 = t2->begin(); it2 != t2->end(); ++it2) {

				Node* v2 = (*it2);

				bitmap C = infos[v2].clade;
				bitmap D = infos[v2].clade_comp;

				
				
				//incompatible iff nb intersections if three or more
				int nbinter = 0;
				if (A.intersects(C))
					++nbinter;
				if (A.intersects(D))
					++nbinter;
				if (B.intersects(C))
					++nbinter;
				if (B.intersects(D))
					++nbinter;

				if (nbinter != 2 && nbinter != 3) {
					ret[v1].push_back(v2);
				}


				//cout << "Comparing v=" << A << "|" << B << "  w=" << C << "|" << D << "      nbi=" << nbinter << endl;

			}
		}

		return ret;
	}

	
};










//that function is chatGPT
map<string, string> parseArguments(int argc, char* argv[]) {
	map<string, string> args;

	for (int i = 1; i < argc; ++i) {
		string arg = argv[i];

		// Check if argument starts with "--" or "-"
		if (arg.rfind("--", 0) == 0 || arg.rfind("-", 0) == 0) {
			// Remove the leading dashes
			string argName = arg.substr(arg.find_first_not_of('-'));

			// Check if the next argument exists and doesn't start with "-"
			if (i + 1 < argc && argv[i + 1][0] != '-') {
				args[argName] = argv[i + 1];
				++i;  // Skip the next argument, as it's the value for the current argument
			}
			else {
				args[argName] = "";  // For flags without a value
			}
		}
	}

	return args;
}




void exec_all_reroots(map<string, string>& args) {
	string infilename = "";
	if (args.count("i"))
		infilename = args["i"];
	else {
		cout << "Please specify an input filename with -i [filename]" << endl;
		return;
	}

	string outfilename = "";
	if (args.count("o"))
		outfilename = args["o"];

	string incontent = Util::GetFileContent(infilename);
	

	if (incontent == "") {
		cout << "Could not find newick string.  Make sure the specified file exists and is non-empty." << endl;
		return;
	}

	Node* root = NewickLex::ParseNewickString(incontent);
	vector<Node*> nodes = root->get_postordered_nodes();

	string outstr = "";
	
	for (Node* v : nodes) {

		if (v->is_root())
			continue;

		Node* w = TreeUtil::subdivide_parent_edge(v);
		TreeUtil::reroot_on_node(w);

		string newick = NewickLex::ToNewickString(w, true, true);
		outstr += newick + "\n";

		TreeUtil::reroot_on_node(root);
		TreeUtil::contract_parent_edge(w);
	}

	if (outfilename == "")
		cout << outstr;
	else {
		Util::WriteFileContent(outfilename, outstr);
	}


}





int main(int argc, char** argv) {

	BipartiteMWIS bwis;

	map<string, string> args = parseArguments(argc, argv);

	/*
	//just some tests
	args["m"] = "rnd";
	args["n"] = "10000";
	//args["o"] = "C:\\Users\\Manuel\\Desktop\\tmp\\trees.txt";

	args["m"] = "all_reroots";
	//args["i"] = "C:\\Users\\Manuel\\Desktop\\tmp\\tree.txt";
	args["i"] = "C:\\Users\\lafm2722\\Desktop\\tmp\\tree.txt";
	*/

	if (args.count("m") && args["m"] == "all_reroots") {
		exec_all_reroots(args);
	}



	if (args.count("m") && args["m"] == "rnd") {
		
		string outfile = "";
		int nbtrees = 2;
		int nbleaves = 10;

		if (args.count("t"))
			nbtrees = Util::ToInt(args["t"]);

		if (args.count("n"))
			nbleaves = Util::ToInt(args["n"]);

		ofstream outfile_stream;
		if (args.count("o")) {
			outfile = args["o"];
			outfile_stream.open(outfile);
		}

		srand(time(NULL));
		
		vector<Node*> trees;

		for (int i = 0; i < nbtrees; ++i) {
			Node* v = new Node();
			TreeUtil::get_random_binary_tree(v, nbleaves);
			
			//contract any edge to unroot
			TreeUtil::contract_parent_edge(v->get_child(0));

			TreeUtil::randomize_branch_lengths(v, 1.0, 10000.0);
			trees.push_back(v);

			string nw = NewickLex::ToNewickString(v, true);
			nw = Util::ReplaceAll(nw, " ", "");




			if (outfile == "")
				cout << nw << endl;
			else
				outfile_stream << nw << endl;
		}




		TreePairInfo tpi(trees[0], trees[1]);
		auto incomp = tpi.get_imcompats();

		cout << "nb incompat=" << incomp.size();
		/*for (auto it = incomp.begin(); it != incomp.end(); ++it) {
			Node* v = (*it).first;
			vector<Node*>& incs = (*it).second;

			for (Node* w : incs) {
				cout << "Incompat" << endl 
					 << "v=" << tpi.infos[v].clade << "|" << tpi.infos[v].clade_comp << endl
					 << "w=" << tpi.infos[w].clade << "|" << tpi.infos[w].clade_comp << endl;
			}
		}*/



		if (outfile != "")
			outfile_stream.close();

		for (Node* v : trees)
			delete v;


	}
	

	return 0;
}



