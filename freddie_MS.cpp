#include "robin_hood.h"
#include "strict_fstream.hpp"
#include "sparse_map.h"
#include "zstr.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>



using namespace std;



typedef __uint128_t kmer;
// typedef uint64_t kmer;
namespace std {
template <>
  struct hash<__uint128_t>
  {
    size_t operator()(const kmer& k) const
    {

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      return (k>>64)+(uint64_t)k;
    }
  };
}

typedef uint8_t color;
// key = kmer, value = pair (color: genomes where the kmer occurs. Max 8 genomes, boolean: has been seen (for feeding venn diagrams))
// typedef tsl::sparse_map<kmer, pair<color,bool>> Map;
typedef robin_hood::unordered_flat_map<kmer, pair<color,bool>> Map;




// severe limitation here, todo: authorize more than 8 colors.
int max_color(8);
int color_number(0);
int k(63);
uint steps(2);
kmer offsetUpdateAnchors = 1;
array<mutex, 1024> nutex;



string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}


void set_color(color& c, int indice) {
	c |= (1 << indice);
}



bool is_set(int indice, color c) {
	c >>= indice;
	return c % 2 == 1;
}


string get_color_code(color c, int n=color_number){
	string res="";
	for(int i(0);i<n;++i){
		// cout<<c%2;
		if (c%2) res+= '1';
		else  res+= '0';
		c>>=1;
	}
	return res;
}


string print_color(color c, int n){
	string res="";
	for(int i(0);i<n;++i){
		cout<<c%2;
		res += c%2;
		c>>=1;
	}
	cout<<" "<<flush;
	return res;
}


// checks if c2 is included in c1. "included" means that all 1's in c1 are also in c2.
// this is not the case on this example as the second bit is 1 in c1 and 0 in c2
// 01010110 c1
// 00110111 c2
bool is_included(color c1, color c2) {
	for (int i(0); i < max_color; ++i) {
		if (c2 % 2 == 0 and c1 % 2 == 1) {
			return false;
		}
		c2 >>= 1;
		c1 >>= 1;
	}
	return true;
}



// optimisable (see eg BankBinary from GATB)
kmer str2num(const string& str) {
	kmer res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		switch (str[i]) {
			case 'A':
				res += 0;
				break;
			case 'C':
				res += 1;
				break;
			case 'G':
				res += 2;
				break;
			default:
				res += 3;
				break;
		}
	}
	return res;
}


/**
 * Reverse complement of a kmer.
 */
kmer rcb(kmer min, uint n) {
	kmer res(0);
	kmer offset(1);
	offset <<= (2 * n - 2);
	for (uint i(0); i < n; ++i) {
		res += (3 - (min % 4)) * offset;
		min >>= 2;
		offset >>= 2;
	}
	return res;
}



kmer nuc2int(char c) {
	switch (c) {
		//~ case 'A': return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
	}
	return 0;
}



kmer nuc2intrc(char c) {
	switch (c) {
		case 'A':
			return 3;
		case 'C':
			return 2;
		case 'G':
			return 1;
			//~ case 'T': return 0;
	}
	return 0;
}


/**
 * Add a nucleotide to an already computed kmer
 */
void updateK(kmer& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchors;
}


/**
 * Add a nucleotide to an already computed reverse complement of a kmer
 */
void updateRCK(kmer& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * (k - 1)));
}


/**
 * hash value for a kmer
 */
kmer hash64shift(kmer key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}




/**
 * Index all canonical kmers from the fasta file (file_name)
 * map then contains for each hashed kmer, the canonical value and its color, updated with this reference_number
 * returns the number of read kmers.
 */
uint64_t load_reference(Map map[], const string file_name, int reference_number, uint i) {
	uint64_t result(0);// nb indexed kmers for this reference.
	if (reference_number >= max_color) {
		cout << "Too much references, max is: " << max_color << endl;
		cout << "I ignore " << file_name << endl;
		return 0;
	}
	zstr::ifstream in(file_name);
	if (not in.good()) {
		cout << "Problem with ref file opening:" << file_name << endl;
		return 0;
	}
	#pragma omp parallel
	while (not in.eof()) {
		string ref, useless;
        #pragma omp critical(file_ref)
		{
			getline(in, useless);   // read a comment, useless
			getline(in, ref);		// read the ACGT sequence
		}
		if (ref.size()>=k) {
			// read all kmers from the ref sequence
			kmer seq(str2num(ref.substr(0, k))), rcSeq(rcb(seq, k)), canon(min(seq, rcSeq));
			uint H(hash64shift(canon));
			uint Hache( H% 16);
            if((H/16)%steps==i){
                nutex[Hache].lock();
                if (map[Hache].count(canon) == 0) {
					// for this kmer, set its indexed value adding the reference_number
					set_color(map[Hache][canon].first, reference_number);
					#pragma omp atomic
					result++;
				}else{
					if(not is_set(reference_number,map[Hache][canon].first)){
						// for this kmer, set its indexed value adding the reference_number
						set_color(map[Hache][canon].first, reference_number);
						#pragma omp atomic
						result++;
					}
				}
				nutex[Hache].unlock();
            }
            
			for (uint64_t j(0); j + k < ref.size(); ++j) {
				updateK(seq, ref[j + k]);
				updateRCK(rcSeq, ref[j + k]);
				canon = (min(seq, rcSeq));
				uint H(hash64shift(canon));
                uint Hache( H% 16);
                if((H/16)%steps==i){
                    nutex[Hache].lock();
                    if (map[Hache].count(canon) == 0) {
						// for this kmer, set its indexed value adding the reference_number
						set_color(map[Hache][canon].first, reference_number);
						#pragma omp atomic
						result++;
					}else{
						if(not is_set(reference_number,map[Hache][canon].first)){
							// for this kmer, set its indexed value adding the reference_number
							set_color(map[Hache][canon].first, reference_number);
							#pragma omp atomic
							result++;
						}
					}
                    nutex[Hache].unlock();
                    
                }
			}
		}
	}
	return result;
}



/**
 * Read the file of file referencing fasta files
 * Returns a vector of the number of read kmers per file
 */
vector<uint64_t> load_reference_file(Map map[], const string& file_name,uint i) {
	vector<uint64_t> result;
	zstr::ifstream in(file_name);
	if (not in.good()) {
		cout << "Problem with ref file opening:" << file_name << endl;
		exit(1);
	}
	string ref_file;
	int reference_number(0);
	while (not in.eof()) {
		getline(in, ref_file);
		if (ref_file.size() > 1) {
			result.push_back(load_reference(map, ref_file, reference_number,i));
            if(i==0){
                color_number++;
            }
            reference_number++;
		}
	}
	return result;
}



/**
 * For each potential color, store the number of corresponding kmers
 */
vector<uint64_t> Venn_evaluation(Map map[], const string& file_name, int size_result,uint i) {
	vector<uint64_t> result(size_result);
	zstr::ifstream in(file_name);
	if (not in.good()) {
		cout << "Problem with ref file opening:" << file_name << endl;
		exit(1);
	}
	#pragma omp parallel
	while (not in.eof()) {
		string ref, useless;
		#pragma omp critical(file_ref)
		{
			getline(in, useless);
			getline(in, ref);
		}
		if (ref.size()>=(uint)k) {
			kmer seq(str2num(ref.substr(0, k)));
			kmer rcSeq(rcb(seq, k));
			kmer canon(min(seq, rcSeq));
            uint H(hash64shift(canon));
			uint Hache( H% 16);
			color c(0);
			bool done=false;
            if((H/16)%steps==i){
                nutex[Hache].lock();
                // enables to store in results once all kmers.
                // optimisable (eg with a set once initialy reading kmers)
                if (map[Hache].count(canon) != 0) { 
                    c = map[Hache][canon].first;
                    done=map[Hache][canon].second;
                    map[Hache][canon].second=true; // validate if this kmer has already been seen
                }
                nutex[Hache].unlock();
                if(not done){
                    #pragma omp atomic
                    ++result[c];
                }
            }

			for (uint64_t j(0); j + k < ref.size(); ++j) {
				updateK(seq, ref[j + k]);
				updateRCK(rcSeq, ref[j + k]);
				canon = (min(seq, rcSeq));
				uint H(hash64shift(canon));
			    uint Hache( H% 16);
                if((H/16)%steps==i){
                    nutex[Hache].lock();
                    if (map[Hache].count(canon) != 0) {
                        c = map[Hache][canon].first;
                        done=map[Hache][canon].second;
                        map[Hache][canon].second=true;
                    }else{
                        c=0;
                        done=false;
                    }
                    nutex[Hache].unlock();
                    if(not done){
                        #pragma omp atomic
                        ++result[c];
                    }
                }
			}
		}
	}
	return result;
}


/**
 * prints the cardinality of each color
 */
void evaluate_completness(const vector<uint64_t>& cardinalities, const vector<uint64_t>& venn, int size_result, const string& venout_file_name = "venn_out.txt") {

	ofstream out(venout_file_name);
	if (not out.good()) {
		cout << "Problem opening file:" << venout_file_name << " for writing the venn resutls"<<endl;
		exit(1);
	}

	vector<uint64_t> counted(cardinalities.size());
	cout<<"Venn:	"<<endl;
	out<<"#Venn:	"<<endl;
	for (uint64_t i(0); i < venn.size(); ++i) {
		string colors = get_color_code(i,size_result);
		cout<<colors<<" "<<intToString(venn[i])<<endl;
		out<<colors<<" "<<venn[i]<<endl;
		uint64_t i_bin(i);
		uint64_t id(0);
		while (i_bin != 0) {
			if (i_bin % 2 == 1) {
				counted[id] += venn[i];
			}
			i_bin >>= 1;
			++id;
		}
	}
	out.close();
	for (uint64_t i(0); i < counted.size(); ++i) {
		cout<<"Kmers  found from file:	" << i << "	" <<intToString(counted[i])  << endl;
		cout<<"Card  of file:	" << i << "	" << cardinalities[i]  << endl;
		cout<<"Completness % for file:	" << i << "	" << (double)100 * counted[i] / cardinalities[i] << endl;
	}
}







int main(int argc, char** argv) {
	if (argc < 4) {
		cout << "[Reference file of file] [query file] [kmer size] (steps number default:2)" << endl;
		exit(0);
	}
	if (argc > 4) {
		steps = (stoi(argv[4]));
	}
	k=(stoi(argv[3]));
	offsetUpdateAnchors <<= (2 * (k));
	auto start = chrono::system_clock::now();
	string inputFILE(argv[1]);
	string inputRef(argv[2]);

	Map map[16];
    vector<uint64_t> venn,cardinalities;
	cout<<"K="+to_string(k)<<endl;
    for(uint i(0);i<steps;++i){
        cout<<"Step "+to_string(i+1)<<endl;
        cout<<"LOAD REFERENCES"<<endl;
        vector<uint64_t> partial_cardinalities(load_reference_file(map, inputFILE,i));
        cout<<"VENN EVALUATION"<<endl;
        vector<uint64_t> partial_venn(Venn_evaluation(map, inputRef,1<<(partial_cardinalities.size()),i));
        if(venn.size()<partial_venn.size()){
            venn.resize(partial_venn.size(),0);
        }
        if(cardinalities.size()<partial_cardinalities.size()){
            cardinalities.resize(partial_cardinalities.size(),0);
        }
        for(uint j(0);j<partial_venn.size();++j){
            venn[j]+=partial_venn[j];
        }
        for(uint j(0);j<partial_cardinalities.size();++j){
           cardinalities[j]+=partial_cardinalities[j];
        }
        for(uint j(0);j<16;++j){
			map[j].clear();
		}
    }
    cout<<cardinalities.size()<<" "<<venn.size()<<endl;
	cout<<"Completness EVALUATION"<<endl;
	evaluate_completness(cardinalities, venn, cardinalities.size());

	auto end                                 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	time_t end_time                          = chrono::system_clock::to_time_t(end);

	cout << "\nFinished computation at " << ctime(&end_time) << "Elapsed time: " << intToString(elapsed_seconds.count()) << "s\n";
}
