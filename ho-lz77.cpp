/*
 * Build co-BWT (using co-lex order) of the text and compute ho-LZ77 pairs
 */

#include <dynamic.hpp>

using namespace dyn;
using namespace std;

void help(){

	cout << "Usage: ho-lz77 input.txt\n";
	exit(0);

}

bool empty_interval(pair<uint64_t, uint64_t> interval){
	return interval.first >= interval.second;
}

/*
 * input: bwt, index on L column, number of steps L
 * output: perform L LF steps from index and return resulting index
 *
 * note: we actually call function FL (not LF) because the underlying BWT
 * is logically suffix-sorted
 */
uint64_t jump_back(wt_bwt & bwt, uint64_t index, uint64_t L){

	for(uint64_t i = 0; i<L; ++i)
		index = bwt.FL(index);

	return index;

}

inline void output_phrase(	wt_bwt & bwt,
							char & c,
							vector<pair<int64_t, uint64_t> > & LZ77k,
							pair<uint64_t,uint64_t> & prev_range,
							pair<uint64_t,uint64_t> & range,
							uint64_t & index,
							string & phrase){

	range = prev_range;

	//closest co-lex positions before and
	//after wt.get_terminator_position()
	//that are followed by the current phrase

	//assign a default high value to previous/next occurrence of the phrase in co-lex order;
	//if they are not initialized later, the high value will not be chosen since
	//we minimize the distance with the current prefix in co-lex order
	int64_t pred_occ = bwt.size()*2;
	int64_t succ_occ = bwt.size()*2;

	//if previous position is valid, compute it
	if(index>0 && index-1 >= range.first && index-1 < range.second){

		pred_occ = jump_back(bwt, index-1, phrase.length());

	}

	//if following position is valid, compute it
	if(index >= range.first && index < range.second){

		succ_occ = jump_back(bwt, index, phrase.length());

	}

	//at least one of the two must be initialized since there must be a previous
	//occurrence of the phrase!
	assert(pred_occ < bwt.size() or succ_occ < bwt.size());

	//co-lex position of current prefix
	int64_t current_prefix_pos = bwt.get_terminator_position();

	//the previous occurrence cannot be the current text prefix!
	assert(pred_occ != current_prefix_pos and succ_occ != current_prefix_pos);

	//absolute distances
	int64_t abs_dist_pred = current_prefix_pos > pred_occ ? current_prefix_pos - pred_occ : pred_occ - current_prefix_pos;
	int64_t abs_dist_succ = current_prefix_pos > succ_occ ? current_prefix_pos - succ_occ : succ_occ - current_prefix_pos;

	int64_t occ = abs_dist_pred < abs_dist_succ ? current_prefix_pos - pred_occ : current_prefix_pos - succ_occ;

	//create new phrase
	LZ77k.push_back({occ,phrase.length()});

	//c remains the previous one that we read, since it is not part of the phrase

	//extend BWT with the characters in the phrase
	for(auto c : phrase)
		bwt.extend(c);

	//erase phrase content
	phrase = string();

	//re-initialize range to full range
	range = bwt.get_full_interval();

	//re-initialize index of current text prefix to (new) position of terminator
	index = bwt.get_terminator_position();

}

int main(int argc,char** argv){

	if(argc<2) help();

	string filePath = argv[1];

	cout << "Input file " << filePath << endl;

	// Build Huffman-encoded BWT

	wt_bwt bwt;
	set<char> alphabet;

	{
		ifstream in(filePath);
		auto F = get_frequencies(in);

		for(auto f : F) if(f.second>0) alphabet.insert(f.first);

		bwt = wt_bwt(F); //Huffman-encoded BWT
	}

	// prepend the alphabet to the text

	string prefix;

	for (auto rit = alphabet.rbegin(); rit != alphabet.rend(); rit++){

		auto c = *rit;
		prefix += c;
		bwt.extend(c);

	}

	// process the text

	vector<pair<int64_t, uint64_t> > LZ77k;// the parse

	{
		ifstream in(filePath);

		char c;
		auto range = bwt.get_full_interval();// interval of current phrase
		pair<uint64_t,uint64_t> prev_range;
		uint64_t index = bwt.get_terminator_position();// co-lex position where current phrase should be if inserted
		string phrase; // current phrase

		if(in.is_open()) {

			if(in.good())
				in.get(c);//get character

			while(in.good()) {

				prev_range = range;
				range = bwt.LF(range, c);

				if(empty_interval(range)){//end of phrase

					output_phrase(bwt,c,LZ77k,prev_range,range,index,phrase);

				}else{

					index = bwt.LF(index,c);
					phrase += c;
					in.get(c);//get next character

				}

			}

		}

		//last phrase has not been output
		if(phrase.length()>0){

			prev_range = range;
			output_phrase(bwt,c,LZ77k,prev_range,range,index,phrase);

		}

		if(!in.eof() && in.fail())
			cout << "Error reading " << filePath << endl;

	}

	cout << "factorization: " << endl;
	for(auto p : LZ77k){

		cout << p.first << ", " << p.second << endl;

	}

	return 0;

}
