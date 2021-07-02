#include <math.h>
#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "sparse_map.h"
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



typedef string kmer;
typedef uint8_t color;

hash<string> hasher;
// key = kmer, value = pair (color: genomes where the kmer occurs. Max 8 genomes, boolean: has been seen (for feeding venn diagrams))
typedef tsl::sparse_map<kmer, pair<color,bool>> Map;



// severe limitation here, todo: authorize more than 8 colors.
uint max_color(8);
uint color_number(0);
uint k(63);
uint nb_bytes(floor((double)k/4));
array<mutex, 1024> nutex;
uint core_number(16);
uint steps(1);



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



void print_char_nuc(uint8_t c) {
    for(uint i=0; i<4;++i){
        switch (c%4) {
			case 0:
				cout<<"A";
				break;
			case 1:
				cout<<"C";
				break;
			case 2:
				cout<<"G";
				break;
			default:
				cout<<"T";
				break;
		}
        c>>=2;
    }
    
}


void print_char_bits(uint8_t c) {
    string toprint;
    for(uint i=0; i<8;++i){
        if(c%2 == 0){
            toprint.push_back('0');
        }else{
            toprint.push_back('1');
        }
        c>>=1;
    }
    reverse(toprint.begin(), toprint.end());
    cout<<toprint;
}


void print_string_bits(const string& str) {
    for(uint i=0; i<str.size();++i){
        print_char_bits(str[i]);
    }
    cout<<endl;
}




string num2str(const kmer& str,const int k) {
    string result;
	for (uint64_t i(0); i < k;++i ) {
		uint8_t byte(str[i/4]);
        byte>>=(2*(3-i%4));
		switch (byte%4) {
			case 0:
				result.push_back('A');
				break;
			case 1:
				result.push_back('C');
				break;
			case 2:
				result.push_back('G');
				break;
			default:
				result.push_back('T');
				break;
		}
	}
	return result;
}



kmer str2num(const string& str) {
    // cout<<str<<endl;
    kmer result;
	unsigned char byte(0);
    uint64_t i;
	for (i=(0); i < str.size(); ) {
		byte <<= 2;
		switch (str[i]) {
			case 'A':
				byte += 0;
				break;
			case 'C':
				byte += 1;
				break;
			case 'G':
				byte += 2;
				break;
			default:
				byte += 3;
				break;
		}
        i++;
        if(i%4==0){
            // print_char_bits(byte);
            // cout<<"push_back"<<endl;
            result.push_back(byte);
            byte = 0;
        }
	}
    if((i%4)!=0){
        byte<<=(2*(4-(i%4)));
        result.push_back(byte);
    }
    if(num2str(result,63)!=str){
        cout<<str<<endl;
        print_string_bits(result);
        cout<<num2str(result,63)<<endl;
        cout<<"probleme"<<endl;
        cin.get();
    }
	return result;
}





void RC_string(string& str){
    reverse(str.begin(),str.end());
    for(int i(0);i<str.size();++i){
        switch (str[i]){
            case 'A':
                str[i] = 'T';
                break;    
            case 'C':
                str[i] = 'G';
                break;
            case 'G':
                str[i] = 'C';
                break;
            default:
                str[i] = 'A';
                break;
        }
    }
    // return str;
}


kmer RC_kmer(const kmer& came, const int k){
    string str(num2str(came,k));
	// string strcopy(str);
	// RC_string(strcopy);
	// RC_string(strcopy);
	// if(str!=strcopy){
	// 	cout<<"RC FAIL"<<endl;
	// 	cin.get();
	// }
    RC_string(str);

    return str2num(str);
}


kmer get_canon(const kmer& str,const int k){
    return min(str, RC_kmer(str,k));    
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
	#pragma omp parallel num_threads(core_number)
	while (not in.eof()) {
        string canon;
		string ref, useless;
    #pragma omp critical(file_ref)
		{
			getline(in, useless);   // read a comment, useless
			getline(in, ref);		// read the ACGT sequence
		}
		if (not ref.empty()) {
			// read all kmers from the ref sequence
			for (uint64_t j(0); j + k < ref.size(); ++j) {
				canon = get_canon(ref.substr(j,k),k);
				uint H(hasher(canon)); 
				uint Hache(H%16);
				if((H/16)%steps==i){
					nutex[Hache].lock();
					// for this kmer, set its indexed value adding the reference_number
					set_color(map[Hache][canon].first, reference_number);
					nutex[Hache].unlock();
					#pragma omp atomic
						result++;
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
			if(i==0) {
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
vector<uint64_t> Venn_evaluation(Map map[], const string& file_name, int size_result, uint i) {
	vector<uint64_t> result(size_result);
	zstr::ifstream in(file_name);
	if (not in.good()) {
		cout << "Problem with ref file opening:" << file_name << endl;
		exit(1);
	}
	#pragma omp parallel num_threads(core_number)
	while (not in.eof()) {
		string ref, useless;
		#pragma omp critical(file_ref)
		{
			getline(in, useless);
			getline(in, ref);
		}
		if (ref.size()>=(uint)k) {

			color c(0);
			bool done=false;
			// enables to store in results once all kmers.
			// optimisable (eg with a set once initialy reading kmers)
            string canon;
			for (uint64_t j(0); j + k < ref.size(); ++j) {
				canon = get_canon(ref.substr(j,k),k);
				uint H(hasher(canon)); 
				uint Hache(H%16);
				if((H/16)%steps==i){
					nutex[Hache].lock();
					if (map[Hache].count(canon) != 0) {
						c = map[Hache][canon].first;
						done=map[Hache][canon].second;
						map[Hache][canon].second=true;
					}else{
						// cout<<"weird"<<endl;
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
		// cout<<"Card  of file:	" << i << "	" << cardinalities[i]  << endl;
		cout<<"Completness % for file:	" << i << "	" << (double)100 * counted[i] / cardinalities[i] << endl;
	}
}





int main(int argc, char** argv) {
	if (argc < 3) {
		cout << "[Reference file of file] [query file] [kmer size] (Steps default:1) (core to use default:16) " << endl;
		exit(0);
	}
	if (argc > 3) {
		k = (stoi(argv[3]));
	}
    if (argc > 4) {
		steps = (stoi(argv[4]));
	}
    if (argc > 5) {
		core_number = (stoi(argv[4]));
	}
	auto start = chrono::system_clock::now();
	string inputFILE(argv[1]);
	string inputRef(argv[2]);

	Map map[16];
	vector<uint64_t> cardinalities,venn;
    for(uint i(0); i < steps;++i){
		cout<<"Step "<<to_string(i)<<endl;
		cout<<"LOAD REFERENCES"<<endl;
		vector<uint64_t> partial_cardinalities(load_reference_file(map, inputFILE,i));
		cout<<"VENN EVALUATION"<<endl;
		vector<uint64_t> partial_venn(Venn_evaluation(map, inputRef,1<<(partial_cardinalities.size()),i));
		if(cardinalities.size()<partial_cardinalities.size()){
			cardinalities.resize(partial_cardinalities.size(),0);
		}
		for(uint j(0);j<partial_cardinalities.size();++j){
			cardinalities[j]+=partial_cardinalities[j];
		}
		if(venn.size()<partial_venn.size()){
			venn.resize(partial_venn.size(),0);
		}
		for(uint j(0);j<partial_venn.size();++j){
			venn[j]+=partial_venn[j];
		}
		for(uint j(0);j<16;++j){
			map[j].clear();
		}
    }
	cout<<"Completness EVALUATION"<<endl;
	evaluate_completness(cardinalities, venn, cardinalities.size());
	auto end                                 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	time_t end_time                          = chrono::system_clock::to_time_t(end);

	cout << "\nFinished computation at " << ctime(&end_time) << "Elapsed time: " << intToString(elapsed_seconds.count()) << "s\n";
}
