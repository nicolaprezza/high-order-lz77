/*
 * Build co-BWT (using co-lex order) of the text and compute ho-LZ77 pairs
 */

#include <dynamic.hpp>
#include <unistd.h>
#include <unordered_map>

using namespace dyn;
using namespace std;

void help(){

	cout << "ho-lz77 [options]" << endl <<
	"Options:" << endl <<
	"-i <arg>   Input file (REQUIRED)" << endl <<
	"-t         Output in triples format (default: pairs)." << endl <<
	"-h         This help." << endl;
	 exit(0);

}
/*
 * number of bits required to write down x>0
 */
uint64_t bit_size(uint64_t x){

	return x==0 ? 1 : 64 - __builtin_clzll(x);

}

/*
 * compute the bit-size of the gamma encoding of x
 */
uint64_t gamma(uint64_t x){

	return 2*bit_size(x) - 1;

}

/*
 * compute the bit-size of the delta encoding of x
 */
uint64_t delta(uint64_t x){

	auto bits = bit_size(x);//bits needed to encode x

	return gamma(bits) + bits -1;

}

/*
 * compute how many bits would the parse take in practice
 */
uint64_t compute_gamma_bit_complexity(vector<pair<int64_t,uint64_t> > & parse){

	uint64_t x = 0;

	for(auto p:parse){

		uint64_t off = p.first<0?-p.first:p.first;//absolute value
		uint64_t len = p.second;

		x++;//1 bit for the sign of off

		x += gamma(off);
		x += gamma(len);

	}

	return x;

}

uint64_t compute_delta_bit_complexity(vector<pair<int64_t,uint64_t> > & parse){

	uint64_t x = 0;

	for(auto p:parse){

		uint64_t off = p.first<0?-p.first:p.first;//absolute value
		uint64_t len = p.second;

		x++;//1 bit for the sign of off

		x += delta(off);
		x += delta(len);

	}

	return x;

}

uint64_t compute_Huff_bit_complexity(vector<pair<int64_t,uint64_t> > & pairs, vector<char> trails){

	uint64_t max_encoded_off = 3200000;//max offset that we encode with Huffman. The others use delta(0).delta(offset) (0 receives the longest code)
	uint64_t max_encoded_len = 100000;//max phrase length that we encode with Huffman. The others use delta(0).delta(len) (0 receives the longest code)

	vector<uint64_t> off_freq(max_encoded_off+1);
	vector<uint64_t> len_freq(max_encoded_len+1);
	vector<uint64_t> char_freq(256);

	uint64_t z = pairs.size();

	for(uint64_t i=0;i<z;++i){

		//absolute value
		pairs[i].first = pairs[i].first < 0 ? -pairs[i].first : pairs[i].first;

		if(pairs[i].first==0 or pairs[i].second==0){
			cout << "err: zero entry in a pair: " << pairs[i].first << " , " << pairs[i].second << endl;
			exit(0);
		}

		if(pairs[i].first <= max_encoded_off) off_freq[uint64_t(pairs[i].first)]++;
		if(pairs[i].second <= max_encoded_len) len_freq[pairs[i].second]++;
		char_freq[uint8_t(trails[i])]++;

	}

	double off_sum = 0;
	double len_sum = 0;
	double char_sum = 0;

	for(auto x : off_freq) off_sum += x;
	for(auto x : len_freq) len_sum += x;
	for(auto x : char_freq) char_sum += x;

	vector<pair<uint64_t,double> > off_f(max_encoded_off+1);
	vector<pair<uint64_t,double> > len_f(max_encoded_len+1);
	vector<pair<uint64_t,double> > char_f(256);

	bool p50=false;
	bool p75=false;
	bool p95=false;


	double F = 0;//cumulative
	for(uint64_t i=0;i<max_encoded_off+1;++i){

		off_f[i].first = i;
		off_f[i].second = double(off_freq[i])/off_sum;

		F += off_f[i].second;

		if(F>=0.5 and not p50){ cout << "50% percentile : " << i << endl; p50=true;}
		if(F>=0.75 and not p75){ cout << "75% percentile : " << i << endl; p75=true;}
		if(F>=0.95 and not p95){ cout << "95% percentile : " << i << endl; p95=true;}

	}

	for(uint64_t i=0;i<max_encoded_len+1;++i){

		len_f[i].first = i;
		len_f[i].second = double(len_freq[i])/len_sum;

	}

	for(uint64_t i=0;i<256;++i){

		char_f[i].first = i;
		char_f[i].second = double(char_freq[i])/char_sum;

	}

	auto off_encoding = alphabet_encoder(off_f);
	auto len_encoding = alphabet_encoder(len_f);
	auto char_encoding = alphabet_encoder(char_f);

	/*
	cout << "Huffman lengths for offsets: " << endl;
	for(uint64_t i = 0;i<=500;++i)
		cout << i << ": " << off_encoding.encode(i).size() << endl;

	cout << "Huffman lengths for lengths: " << endl;
	for(uint64_t i = 0;i<=500;++i)
		cout << i << ": " << len_encoding.encode(i).size() << endl;

	cout << "Huffman lengths for trail chars: " << endl;
	for(uint64_t c = 0;c<256;++c)
		cout << char(c) << ": " << char_encoding.encode(uint8_t(c)).size() << endl;*/


	//uint64_t x = (max_encoded_off+2 + max_encoded_len + 2 + 256)*64; //the model
	uint64_t x = 0;

	for(uint64_t i=0;i<z;++i){

		x++; // sign of offset

		if(pairs[i].first <= max_encoded_off)
			x += off_encoding.encode(uint64_t(pairs[i].first)).size();
		else
			x += off_encoding.encode(0).size() + delta(uint64_t(pairs[i].first));

		if(pairs[i].second <= max_encoded_len)
			x += len_encoding.encode(pairs[i].second).size();
		else
			x += len_encoding.encode(0).size() + delta(pairs[i].second);

		x += char_encoding.encode(uint8_t(trails[i])).size();


	}

	return x;

}


/*
 * entropy of a vector
 */
double entropy(vector<uint64_t> & V){

	unordered_map<uint64_t, double> freq;
	double n = V.size();

	for(auto x : V)	freq[x]=0;
	for(auto x : V)	freq[x]++;

	double H = 0;

	for(auto p : freq){

		auto f = p.second/n;

		H -= f*log2(f);

	}

	return H;

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

	//extend BWT with the characters in the phrase
	for(auto c : phrase)
		bwt.extend(uint8_t(c));

	//erase phrase content
	phrase = string();

	//re-initialize range to full range
	range = bwt.get_full_interval();

	//re-initialize index of current text prefix to (new) position of terminator
	index = bwt.get_terminator_position();

}

inline void output_phrase(	wt_bwt & bwt,
							vector<pair<int64_t, uint64_t> > & LZ77k,
							vector<char> & trail_chars,
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

		pred_occ = jump_back(bwt, index-1, phrase.length()-1);

	}

	//if following position is valid, compute it
	if(index >= range.first && index < range.second){

		succ_occ = jump_back(bwt, index, phrase.length()-1);

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
	LZ77k.push_back({occ,phrase.length()-1});
	trail_chars.push_back(phrase[phrase.length()-1]);

	//c remains the previous one that we read, since it is not part of the phrase

	//extend BWT with the characters in the phrase
	for(auto c : phrase)
		bwt.extend(uint8_t(c));

	//erase phrase content
	phrase = string();

	//re-initialize range to full range
	range = bwt.get_full_interval();

	//re-initialize index of current text prefix to (new) position of terminator
	index = bwt.get_terminator_position();

}

void run_pairs(string filePath){

	cout << "Computing the parse in format (occ,len)" << endl;

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
		bwt.extend(uint8_t(c));

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

			uint64_t read_char = 1;

			while(in.good()) {

				prev_range = range;
				range = bwt.LF(range, uint8_t(c));

				if(empty_interval(range)){//end of phrase

					output_phrase(bwt,LZ77k,prev_range,range,index,phrase);

				}else{

					index = bwt.LF(index,uint8_t(c));
					phrase += c;
					in.get(c);//get next character
					read_char++;

					if(read_char%100000 == 0){

						cout << "read " << read_char << " characters." << endl;

					}

				}

			}

		}

		//last phrase has not been output
		if(phrase.length()>0){

			prev_range = range;
			output_phrase(bwt,LZ77k,prev_range,range,index,phrase);

		}

		if(!in.eof() && in.fail())
			cout << "Error reading " << filePath << endl;

	}

	/*cout << "factorization: " << endl;
	for(auto p : LZ77k){

		cout << p.first << ", " << p.second << endl;

	}*/

	auto N = bwt.size() -1;//file length


	uint64_t positive = 0;//positive offsets

	for(auto p : LZ77k){

		positive += p.first>0;

	}

	cout << "positive offsets: " << positive << endl;
	cout << "negative offsets: " << LZ77k.size()-positive << endl;

	int bucket_size = 1;

	auto buckets = vector<uint64_t>(bwt.size()/bucket_size + 1);

	for(auto p : LZ77k){

		buckets[(p.first<0?-p.first:p.first)/bucket_size]++;

	}

	for(int i=0;i<1000;++i){

		//cout << "[" << i*bucket_size << "," << (i+1)*bucket_size << ") : " << buckets[i] << endl;
		cout << i << "\t" << buckets[i] << endl;

	}

	cout << "number of phrases = " << LZ77k.size() << endl;
	cout << "gamma complexity of the output: " << compute_gamma_bit_complexity(LZ77k)/8+1 << " Bytes, " << double(compute_gamma_bit_complexity(LZ77k))/double(N) << " bit/symbol" << endl;
	cout << "delta complexity of the output: " << compute_delta_bit_complexity(LZ77k)/8+1 << " Bytes, " << double(compute_delta_bit_complexity(LZ77k))/double(N) << " bit/symbol" << endl;

}


void run_triples(string filePath){

	cout << "Computing the parse in format (occ,len,char)" << endl;

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
		bwt.extend(uint8_t(c));

	}

	// the parse
	vector<pair<int64_t, uint64_t> > LZ77k;
	vector<char> trail_chars;

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

			uint64_t read_char = 1;

			while(in.good()) {

				prev_range = range;
				range = bwt.LF(range, uint8_t(c));

				phrase += c;

				if(empty_interval(range)){//end of phrase

					output_phrase(bwt,LZ77k,trail_chars,prev_range,range,index,phrase);

				}else{

					index = bwt.LF(index,uint8_t(c));

				}

				in.get(c);//get next character
				read_char++;

				if(read_char%100000 == 0){

					cout << "read " << read_char << " characters." << endl;

				}

			}

		}

		//last phrase has not been output
		if(phrase.length()>0){

			output_phrase(bwt,LZ77k,trail_chars,prev_range,range,index,phrase);

		}

		if(!in.eof() && in.fail())
			cout << "Error reading " << filePath << endl;

	}

	/*cout << "factorization: " << endl;
	for(auto p : LZ77k){

		cout << p.first << ", " << p.second << endl;

	}*/

	auto N = bwt.size() -1;//file length


	uint64_t positive = 0;//positive offsets

	for(auto p : LZ77k){

		positive += p.first>0;

	}

	cout << "positive offsets: " << positive << endl;
	cout << "negative offsets: " << LZ77k.size()-positive << endl;

	int bucket_size = 1;

	auto buckets = vector<uint64_t>(bwt.size()/bucket_size + 1);

	for(auto p : LZ77k){

		buckets[(p.first<0?-p.first:p.first)/bucket_size]++;

	}

	/*for(int i=0;i<1000;++i){

		//cout << "[" << i*bucket_size << "," << (i+1)*bucket_size << ") : " << buckets[i] << endl;
		cout << i << "\t" << buckets[i] << endl;

	}*/

	uint64_t gamma_trail = 0;
	uint64_t delta_trail = 0;

	for(auto c : trail_chars){

		gamma_trail += gamma(uint64_t(uint8_t(c))+1);
		delta_trail += delta(uint64_t(uint8_t(c))+1);

	}

	auto G = gamma_trail+compute_gamma_bit_complexity(LZ77k);
	auto D = delta_trail+compute_delta_bit_complexity(LZ77k);
	auto H = compute_Huff_bit_complexity(LZ77k,trail_chars);

	vector<uint64_t> abs_off;
	for(auto x : LZ77k) abs_off.push_back(x.first<0?-x.first:x.first);

	uint64_t sum_log = 0;

	for(auto x:abs_off) sum_log += bit_size(x)+1;

	cout << "number of phrases = " << LZ77k.size() << endl;
	cout << "Entropy of the offsets: " << (entropy(abs_off)+1) << endl;
	cout << "Sum of logs of the offsets: " << sum_log << endl;
	cout << "gamma complexity of the output: " << (G/8)+1 << " Bytes, " << double(G)/double(N) << " bit/symbol" << endl;
	cout << "delta complexity of the output: " << (D/8)+1 << " Bytes, " << double(D)/double(N) << " bit/symbol" << endl;
	cout << "Huffman complexity of the output: " << (H/8)+1 << " Bytes, " << double(H)/double(N) << " bit/symbol" << endl;


}



int main(int argc,char** argv){

	if(argc<2) help();

	bool triples = false;
	string filePath;

	int opt;
	while ((opt = getopt(argc, argv, "i:th")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 't':
				triples = true;
			break;
			case 'i':
				filePath = string(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	if(filePath.length()==0) help();

	cout << "Input file " << filePath << endl;

	if(triples){

		run_triples(filePath);

	}else{

		run_pairs(filePath);

	}

	return 0;

}
