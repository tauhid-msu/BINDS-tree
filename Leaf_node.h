#include "logClass.h"
#include "Node.h"
#include "Leaf_entry.h"
#include <algorithm>
#include <iostream>
const int LEAF_NODE_OVERHEAD = 4; // int count

//extern const int DIR_NODE_OVERHEAD;

extern int debug_boxQ_leaf_accessed;
extern int debug_boxQ_leaf_hit;
extern int debug_boxQ_leaf_hit_peak;
extern logClass logO;

class Leaf_node:public Node {
public:

    //Error_code retrieve_by_box_query(const unsigned char box_query_data_DMBR[DMBR_SIZE], Leaf_entry* query_results, int& number_of_query_results);
    //Error_code retrieve_by_box_query_use_link(const unsigned char box_query_data_DMBR[DMBR_SIZE], Leaf_entry* query_results, int& number_of_query_results,int& number_of_io);
    Error_code retrieve_by_box_query_use_link_v2(const unsigned char box_query_data_DMBR[DMBR_SIZE], Leaf_entry* query_results, int& number_of_query_results,int& number_of_io);
    bool is_within_box(const unsigned char box_query_data_DMBR[DMBR_SIZE], unsigned char* DMBR);
    Leaf_node(int* alphabet_sizes);
    Error_code retrieve(Leaf_entry & query_data);
    Error_code retrieve_by_hamming_dist(const Leaf_entry & query_data, int range, Leaf_entry* query_results, int& number_of_query_results,int &number_of_io);
    Error_code insert_new_data(Leaf_entry &new_data, double leaf_min_util, int *alphabet_sizes, unsigned char* cur_DMBR, Leaf_node* & new_leaf_node, unsigned char* new_DMBR);
    Error_code insert_new_data_sprt(Leaf_entry &new_data, double leaf_min_util, int *alphabet_sizes, unsigned char* cur_DMBR, Leaf_node* & new_leaf_node,
                                    unsigned char* new_DMBR, int& best_dim);
    void set_node_count(int new_count);
    void set_node_entry(int index, Leaf_entry new_entry);
    int get_node_count();
    void get_node_entry(int index, Leaf_entry & cur_entry);
    void get_DMBR(unsigned char* cur_DMBR);
    void get_single_DMBR(unsigned char* cur_DMBR,int index);

    //virtual void read_node(fstream& ND_file, unsigned int block_number);
    //virtual void write_node(fstream& ND_file, unsigned int block_number);
     void read_node(fstream& ND_file, unsigned int block_number);
     void write_node(fstream& ND_file, unsigned int block_number);

    //use link
    static const int LEAF_NODE_SIZE =  (DISK_BLOCK_SIZE - LEAF_NODE_OVERHEAD) / sizeof(Leaf_entry);

    void logDMBR( );
    void log2file( );
    void print_OneEntry_OnOneDscDim(const Leaf_entry*  const oneEntry,int entryIndex, int dimNum);
    void log_OneEntry_OnOneDscDim(const Leaf_entry*  const oneEntry,int entryIndex, int dimNum);

private:
    void split_algorithm_TO_DEATH(Leaf_entry * split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
    void split_TO_DEATH_RANDOM_UTIL(Leaf_entry * split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
    void split_TO_DEATH_RANDOM_UTIL_v3(Leaf_entry * split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
    int split_TO_DEATH_RANDOM_UTIL_v4(Leaf_entry * split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);

    void split_algorithm(Leaf_entry * split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
    int split_algorithm_sprt(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1,
        unsigned char* DMBR2);
    int split_algorithm_after_knapsack(Leaf_entry * split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);

    void create_DMBR(Leaf_entry * leaf_entries, int num_of_leaf_entries, unsigned char* new_DMBR); // Create a DMBR from a set of leaf_entries
    void create_DMBR(Leaf_entry* leaf_entries, int* indices_of_leaf_entries_used, int num_of_leaf_entries_used, unsigned char* new_DMBR); // Create a DMBR from a subset of leaf_entries
    void sort_entries_by_size(int sort_dim, Leaf_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list);
    int sort_entries_by_size_RANDOM_UTIL(int sort_dim, Leaf_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list);
    int sort_entries_by_size_RANDOM_UTIL_v2(int sort_dim, Leaf_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list,double leaf_min_util);
    int sort_entries_by_size_RANDOM_UTIL_v3(int sort_dim, Leaf_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list,double leaf_min_util,int & left_letters);
    void printDMBR(unsigned char* DMBR);
    
    void recursiveUtilCal(int& alphabetSize, int minUtilReq, vector<unsigned char>& varAlphabetFreq, int index, bool* minSplitChars, bool* curSplitChars, int& curMinUtil, int newMinUtil);
    int greedyMinBalanceSplit(int& alphabetSize, unsigned char* alphabet2, vector<unsigned char>& varAlphabetFreq, bool* minSplitChars, double leaf_min_util,
        int curMaxAlphabetSize);
    void NDT_sort_subset(int count, bool asc, vector<unsigned int>& splitDims, vector<unsigned char>& splitDimsCount);
    void calEntriesCharsHistogram(int sort_dim, Leaf_entry* sort_entries, int entriesCount, int& alphabetSize, unsigned char* alphabet2, 
        int* letterFreq, vector<vector<int> >& letter_entry_set, int curMaxAlphabetSize);
    //int sort_entries_by_size_RANDOM_UTIL_enforceUtil(int sort_dim, Leaf_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list,double leaf_min_util);
	//int getCompressedDiskSize(unsigned char* allEntries_DMBR,int numEntries, int *alphabet_sizes);

    // Data members
    int count; // curent number of entries in the node
    Leaf_entry entries[LEAF_NODE_SIZE];

    // used in sort entries by size
    struct letter_entry{
        int freq; // number of entries whose letter corresponds to the letter_entry
        int entry_set[LEAF_NODE_SIZE + 1]; // entries (indices) whose letter corresponds to the letter_entry
    };
};

Leaf_node::Leaf_node(int* alphabet_sizes):Node(alphabet_sizes)
{
   count = 0;
}

void moveToLine(fstream file, int lineNum)
{
}


void update_record(int record_id)
{
fstream readid_file, typeid_file;
const char* readid_filename = (globalRecordFilename+".readid").c_str();
const char* typeid_filename = (globalRecordFilename+".typeid").c_str();

readid_file.open(readid_filename,fstream::in | fstream::out);
typeid_file.open(typeid_filename,fstream::in | fstream::out);

if(readid_file.fail())
{
    cout<<"can't open file .readid"<<endl;
    exit(1);
}
if(typeid_file.fail())
{
    cout<<"can't open file .typeid"<<endl;
    exit(1);
}
unsigned char type_array[256]={'\0'};
int type_num;
typeid_file.seekg(record_id * sizeof(type_array), ios::beg);
typeid_file.read((char*)type_array, sizeof(type_array));

type_num = type_array[0];
int i;
for(i=1; i<=type_num; i++)
{
    if(typeid_global == type_array[i])
	break;
}
if(i<=type_num)
    return;

type_num++;
type_array[type_num] = typeid_global;
type_array[0] = type_num;

/*
int min =0;

for(int i=1; i<type_num+1; i++)
{
    min = i;
    for(int j=i; j<type_num+1; j++)
    {
	if(type_array[j]<type_array[min])
	{
	    min = j;
	}
    }
    if(min != i)
    {
	int tmp = type_array[i];
	type_array[i] = type_array[min];
	type_array[min] = tmp;
    }
    
}
*/
typeid_file.seekp(record_id * sizeof(type_array), ios::beg);

typeid_file.write((const char*)type_array, sizeof(type_array));

readid_file.close();
typeid_file.close();

}

void Leaf_node::printDMBR(unsigned char* DMBR){
  for(int i=0; i<DIM; i++){
    cout<<setw(4)<<( static_cast<unsigned int>( DMBR[ i ]) );
  }
  cout<<endl;
}


void update_record_array(int record_id)
{

int type_num;
type_num = record_type[record_id][0];
int i;
for(i=1; i<=type_num; i++)
{
    if(typeid_global == record_type[record_id][i])
	break;
}
if(i<=type_num)
    return;

type_num++;
record_type[record_id][type_num] = typeid_global;
record_type[record_id][0] = type_num;

}


Error_code Leaf_node::retrieve(Leaf_entry& query_data)
{
    int i, j;
    for(i = 0; i < count; i++)
	{
        for(j = 0; j < DIM; j++)
            if(entries[i].key[j] != query_data.key[j]) 
                break;
        if(j == DIM)
        {
            //query_data = entries[i];           
	    entries[i].record_count++;

update_record(entries[i].record);	    

//array record_type for test
update_record_array(entries[i].record);
//array record_type for test


	    //make a link to record
	    //entries[i].record = query_data.record;

static int vf=0;
vf+=entries[i].record_count;
//cout<<vf<<endl;
//            cout<<entries[i].record<<endl;

//write to file


	    break;
        }
    }
    if(i == count) return not_present;
    else return success;
}

bool Leaf_node::is_within_box(const unsigned char box_query_data_DMBR[DMBR_SIZE], unsigned char* DMBR)
{
    int j,i;
    bool matchOnThisDim;
    for(j = 0; j < DMBR_SIZE; j=j+BYTES_PER_DIM_IN_DMBR)
    {
        matchOnThisDim=false;

        for(i=0;i<BYTES_PER_DIM_IN_DMBR;i++)
            if((DMBR[j+i] & box_query_data_DMBR[j+i])!= 0)
                matchOnThisDim = true;

        if (matchOnThisDim == false)
            return false;
    }
    return true;
}

void output_record(int record_id)
{
    record_id--;
    fstream typeid_file, query_result_file;
    const char* typeid_filename = (globalRecordFilename+".typeid").c_str();
    const char* query_result_filename = (globalBQFilename+ ".result").c_str();
    typeid_file.open(typeid_filename, fstream::in | fstream::out);
    query_result_file.open(query_result_filename, fstream::in | fstream::out| fstream::app);

    if(typeid_file.fail())
    {
	cout<<"can't open file "<<typeid_filename<<endl;
    }
    if(query_result_file.fail())
    {
	cout<<"can't open file "<<query_result_filename<<endl;
    }

    unsigned char type_array[256]={'\0'};
    int type_num;
    typeid_file.seekg(record_id * sizeof(type_array));
    typeid_file.read((char*)type_array, sizeof(type_array));

    type_num = type_array[0];
    type_num++;
    int i;
    for(i=0; i<type_num; i++)
	query_result_file << type_array[i]<<" ";
    query_result_file << endl;

    typeid_file.close();
    query_result_file.close();

}

//dec13
//origin from retrieve_by_box_query_use_link(...)
//use different ways to increase i/o number
Error_code Leaf_node::retrieve_by_box_query_use_link_v2(
    const unsigned char box_query_data_DMBR[DMBR_SIZE],  
    Leaf_entry* query_results, 
    int& number_of_query_results,
    int& number_of_io){
  
    int i;
    int orig_count = number_of_query_results;
    int debug_temp_hit_num=0;

    int totalEntrieNum=0;

    debug_boxQ_leaf_accessed++;
    
    //cout<<"EnCount:"<<count<<endl;
    
    for(i = 0; i < count; i++){

        totalEntrieNum+=entries[i].record_count;

        unsigned char tmp_DMBR[DMBR_SIZE];

        get_single_DMBR(tmp_DMBR,i);

        if(is_within_box(box_query_data_DMBR, tmp_DMBR))
        {
            if(number_of_query_results < QUERY_RESULTS_BUFFER_SIZE)
                query_results[number_of_query_results] = entries[i];

            number_of_query_results++;
            debug_boxQ_leaf_hit++;
            debug_temp_hit_num++;
        }
    }
    int numberOfEntriesInALinkNode = (DISK_BLOCK_SIZE - sizeof(int) )/ (DIM + sizeof(float)*CNTDIM+sizeof(int)/*pointer to database*/);
    //number_of_io += (totalEntrieNum - count)/ numberOfEntriesInALinkNode;
    if(((totalEntrieNum - count) %numberOfEntriesInALinkNode)!=0)
        //number_of_io++;
    if(debug_temp_hit_num>debug_boxQ_leaf_hit_peak)
        debug_boxQ_leaf_hit_peak=debug_temp_hit_num;

    if(orig_count == number_of_query_results) 
        return not_present;
    else
        return success;
}


Error_code Leaf_node::retrieve_by_hamming_dist(const Leaf_entry& query_data, int range, Leaf_entry * query_results, int& number_of_query_results,int& number_of_io)
{
    int i, j;
    int orig_count = number_of_query_results;
    int totalEntrieNum=0;

    int cur_dist;
    for(i = 0; i < count; i++)
    {
        totalEntrieNum+=entries[i].record_count;
        cur_dist = 0;
        for(j = 0; j < DIM && cur_dist <= range; j++)
            if(entries[i].key[j] != query_data.key[j])cur_dist++;
        if(cur_dist <= range)
        { // current vector within range
            if(number_of_query_results < QUERY_RESULTS_BUFFER_SIZE)
                query_results[number_of_query_results] = entries[i];
            number_of_query_results++;
        }
    }
    //if(USE_LINK_FOR_RANGE_QUERY==1)
    //{
    //        int numberOfEntriesInALinkNode = (DISK_BLOCK_SIZE - sizeof(int) )/ (DIM + sizeof(float)*CNTDIM+sizeof(int)/*pointer to database*/);
    //        number_of_io += (totalEntrieNum - count)/ numberOfEntriesInALinkNode;
    //        if(((totalEntrieNum - count) %numberOfEntriesInALinkNode)!=0)
    //            number_of_io++;
    //}
    if(orig_count == number_of_query_results)
        return not_present;
    else
        return success;
}

std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

Error_code Leaf_node::insert_new_data_sprt(Leaf_entry &new_data, double leaf_min_util, int *alphabet_sizes, unsigned char* cur_DMBR, Leaf_node* & new_leaf_node,
unsigned char* new_DMBR, int& best_dim){
  
  /*
  //Find the last line of record file, record is the line number
  fstream record_file, readid_file, typeid_file;
  const char* record_filename = (globalRecordFilename).c_str();
  const char* readid_filename = (globalRecordFilename+".readid").c_str();
  const char* typeid_filename = (globalRecordFilename+".typeid").c_str();
  
  record_file.open(record_filename,fstream::in | fstream::out);
  readid_file.open(readid_filename,fstream::in | fstream::out | fstream::app);
  typeid_file.open(typeid_filename,fstream::in | fstream::out | fstream::app);
  
  if(record_file.fail())
  {
    cout<<"can't open file .record"<<endl;
    exit(1);
  }
  if(readid_file.fail())
  {
    cout<<"can't open file .readid"<<endl;
    exit(1);
  }
  if(typeid_file.fail())
  {
    cout<<"can't open file .typeid"<<endl;
    exit(1);
  }
  
  static int record_count=-1;
  
  record_count++;
  record_file << record_count <<endl;
  
  unsigned char type_array[256]={'\0'};
  type_array[0] = 1;
  type_array[1] = typeid_global;
  //write (read id, type id) to record file
  typeid_file.write((const char*)type_array, sizeof(type_array));
  //typeid_file << "1 " << typeid_global << endl;
  
  new_data.record = record_count;
  record_file.close();
  readid_file.close();
  typeid_file.close();
  
  //array record_type[][] for test
  record_type[record_count][0]=1;
  record_type[record_count][1]=typeid_global;
  //
  
  */
  
  if(count < LEAF_NODE_SIZE)
  { // room avaialbe
    entries[count] = new_data;
    count++;
    // create a new DMBR
    create_DMBR(entries, count, cur_DMBR);
    return success;
  }
  else
  { // overflows, split
    int i;
    Leaf_entry tmp_entries[LEAF_NODE_SIZE+1];
    for(i = 0; i < LEAF_NODE_SIZE; i++)
      tmp_entries[i] = entries[i];
    tmp_entries[LEAF_NODE_SIZE] = new_data;
    int entry_group1[LEAF_NODE_SIZE], count1, entry_group2[LEAF_NODE_SIZE], count2;
    switch (nodeSplitType)
    {
    case ORIGINAL:
      split_algorithm(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      break;
    case SPACE_PARTITION:
      //best_dim = split_algorithm_sprt(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      best_dim = split_TO_DEATH_RANDOM_UTIL_v4(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      break;
    case TO_DEATH:
      split_algorithm_TO_DEATH(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      break;
      //case TO_DEATH_MAX_UTIL:
      //    split_TO_DEATH_MAX_UTIL(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      //    break;
    case TO_DEATH_RANDOM_UTIL:
      //split_TO_DEATH_RANDOM_UTIL(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      if(try_all_dim_with_minUtil==0)
        split_TO_DEATH_RANDOM_UTIL_v3(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      else
        split_TO_DEATH_RANDOM_UTIL_v4(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
      
      break;
    default:
      cout<<"nodeSplitType not defined"<<endl;
    exit(0);
    }
    assert((count1 + count2) == LEAF_NODE_SIZE+1);
    
    // create a new leaf node for group1
    new_leaf_node = new Leaf_node(alphabet_sizes);
    new_leaf_node->set_node_count(count1);
    for(i = 0; i < count1; i++)
      new_leaf_node->set_node_entry(i, tmp_entries[entry_group1[i]]); 
    
    // modify the split (current) node based on group2
    count = count2;
    for(i = 0; i < count2; i++)
      entries[i] = tmp_entries[entry_group2[i]]; // insert entries knowing it will not overflow
    return overflow;
  }
}



Error_code Leaf_node::insert_new_data(Leaf_entry &new_data, double leaf_min_util, int* alphabet_sizes, unsigned char* cur_DMBR, Leaf_node*& new_leaf_node,
unsigned char* new_DMBR){


/*
//Find the last line of record file, record is the line number
fstream record_file, readid_file, typeid_file;
const char* record_filename = (globalRecordFilename).c_str();
const char* readid_filename = (globalRecordFilename+".readid").c_str();
const char* typeid_filename = (globalRecordFilename+".typeid").c_str();

record_file.open(record_filename,fstream::in | fstream::out);
readid_file.open(readid_filename,fstream::in | fstream::out | fstream::app);
typeid_file.open(typeid_filename,fstream::in | fstream::out | fstream::app);

if(record_file.fail())
{
    cout<<"can't open file .record"<<endl;
    exit(1);
}
if(readid_file.fail())
{
    cout<<"can't open file .readid"<<endl;
    exit(1);
}
if(typeid_file.fail())
{
    cout<<"can't open file .typeid"<<endl;
    exit(1);
}

static int record_count=-1;

record_count++;
record_file << record_count <<endl;

unsigned char type_array[256]={'\0'};
type_array[0] = 1;
type_array[1] = typeid_global;
//write (read id, type id) to record file
typeid_file.write((const char*)type_array, sizeof(type_array));
//typeid_file << "1 " << typeid_global << endl;

new_data.record = record_count;
record_file.close();
readid_file.close();
typeid_file.close();

//array record_type[][] for test
record_type[record_count][0]=1;
record_type[record_count][1]=typeid_global;

*/

//
	
    if(count < LEAF_NODE_SIZE)
    { // room avaialbe
        entries[count] = new_data;
        count++;
        // create a new DMBR
        create_DMBR(entries, count, cur_DMBR);
        return success;
    }
    else
    { // overflows, split
        int i;
        Leaf_entry tmp_entries[LEAF_NODE_SIZE+1];
        for(i = 0; i < LEAF_NODE_SIZE; i++)
            tmp_entries[i] = entries[i];
        tmp_entries[LEAF_NODE_SIZE] = new_data;
        int entry_group1[LEAF_NODE_SIZE], count1, entry_group2[LEAF_NODE_SIZE], count2;
        switch (nodeSplitType)
        {
            case ORIGINAL:
                split_algorithm(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
                break;
            case TO_DEATH:
                split_algorithm_TO_DEATH(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
                break;
            //case TO_DEATH_MAX_UTIL:
            //    split_TO_DEATH_MAX_UTIL(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
            //    break;
            case TO_DEATH_RANDOM_UTIL:
                //split_TO_DEATH_RANDOM_UTIL(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
			if(try_all_dim_with_minUtil==0)
                split_TO_DEATH_RANDOM_UTIL_v3(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.
			else
                split_TO_DEATH_RANDOM_UTIL_v4(tmp_entries, leaf_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR.

                break;
            default:
            cout<<"nodeSplitType not defined"<<endl;
            exit(0);
        }
        assert((count1 + count2) == LEAF_NODE_SIZE+1);

        // create a new leaf node for group1
        new_leaf_node = new Leaf_node(alphabet_sizes);
        new_leaf_node->set_node_count(count1);
        for(i = 0; i < count1; i++)
            new_leaf_node->set_node_entry(i, tmp_entries[entry_group1[i]]); 

        // modify the split (current) node based on group2
        count = count2;
        for(i = 0; i < count2; i++)
            entries[i] = tmp_entries[entry_group2[i]]; // insert entries knowing it will not overflow
        return overflow;
    }
}

void Leaf_node::set_node_count(int new_count)
{
    assert(count >= 0);
    count = new_count;
}

void Leaf_node::set_node_entry(int index, Leaf_entry new_entry)
{
    assert(index >= 0 && index < count);
    entries[index] = new_entry;
}

int Leaf_node::get_node_count()
{
    return count;
}

void Leaf_node::get_node_entry(int index, Leaf_entry& cur_entry){
   assert(index >= 0 && index < count);
   cur_entry = entries[index];
}

void Leaf_node::get_single_DMBR(unsigned char* cur_DMBR,int index)
{
    int i, j;
    int byte_no, bit_no;

    for(i = 0; i < DMBR_SIZE; i++)
        cur_DMBR[i] = 0;

    for(j = 0; j < DIM; j++)
    {
        byte_no = this->DMBR_byte_lut[j][entries[index].key[j]];
        bit_no = this->DMBR_bit_lut[j][entries[index].key[j]];
        cur_DMBR[byte_no] |= MASKS[bit_no];
    }
}

void Leaf_node::get_DMBR(unsigned char* cur_DMBR)
{
    int i, j;
    int byte_no, bit_no;

    for(i = 0; i < DMBR_SIZE; i++)
        cur_DMBR[i] = 0;
    for(i = 0; i < count; i++)
    {
        for(j = 0; j < DIM; j++)
        {
            byte_no = this->DMBR_byte_lut[j][entries[i].key[j]];
            bit_no = this->DMBR_bit_lut[j][entries[i].key[j]];
            cur_DMBR[byte_no] |= MASKS[bit_no];
        }
    }
}

void Leaf_node::read_node(fstream& ND_file, unsigned int block_number)
{
    ND_file.seekg(block_number*DISK_BLOCK_SIZE);

    ND_tree_record_count record_count_array[LEAF_NODE_SIZE];
    ND_tree_record record_array[LEAF_NODE_SIZE];
    unsigned char key_array[LEAF_NODE_SIZE][DIM];
    unsigned char cntkey_array[LEAF_NODE_SIZE][CNTDIM];

    ND_file.read((char*)(&count), sizeof(int));

    ND_file.read((char*)record_count_array, sizeof(record_count_array));
    ND_file.read((char*)record_array, sizeof(record_array));
    ND_file.read((char*)key_array, sizeof(key_array));
    ND_file.read((char*)cntkey_array, sizeof(cntkey_array));
    for(int i = 0; i < count; i++)
    {
        entries[i].record_count = record_count_array[i];
	entries[i].record = record_array[i];
        for(int j =0; j < DIM; j++)
            entries[i].key[j] = key_array[i][j];
        for(int j =0; j < CNTDIM; j++)
            entries[i].cntkey[j] = cntkey_array[i][j];
    }
}



void Leaf_node::write_node(fstream& ND_file, unsigned int block_number)
{
    ND_file.seekp(block_number * DISK_BLOCK_SIZE);
    ND_file.write(reinterpret_cast<const char*>(&count), sizeof(int));
    ND_tree_record_count record_count_array[LEAF_NODE_SIZE];
    ND_tree_record record_array[LEAF_NODE_SIZE];
    unsigned char key_array[LEAF_NODE_SIZE][DIM];
    float cntkey_array[LEAF_NODE_SIZE][CNTDIM];

    for(int i = 0; i < count; i++)
    {
        record_count_array[i] = entries[i].record_count;
	record_array[i] = entries[i].record;
	for(int j =0; j < DIM; j++)
            key_array[i][j] = entries[i].key[j];
        for(int j =0; j < CNTDIM; j++)
            cntkey_array[i][j] = entries[i].cntkey[j];
    }

    ND_file.write((const char*)record_count_array, sizeof(record_count_array));
    ND_file.write((const char*)record_array, sizeof(record_array));
    ND_file.write((const char*)key_array, sizeof(key_array));
    ND_file.write((const char*)cntkey_array, sizeof(cntkey_array));
}

//sort all dimensions
//try to split on all dimensions with more than 1 letter
int Leaf_node::split_TO_DEATH_RANDOM_UTIL_v4(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
    int i, j, k;

    int num_of_entries = LEAF_NODE_SIZE + 1;
    int sorted_entry_list[LEAF_NODE_SIZE + 1];

    //int split_start = static_cast<int>(floor(LEAF_NODE_SIZE * leaf_min_util)) - 1; // the index in the entries where possible split point starts
    //if(split_start < 0)
    //    split_start = 0;
    //int split_end = static_cast<int>(ceil(LEAF_NODE_SIZE * (1 - leaf_min_util))); // the index in the entries where possible split point ends
    //if(split_end >= num_of_entries - 1)
    //    split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
    //if(split_end < split_start)
    //    split_end = split_start;
    count1 = 0;
    count2 = 0;

    // sort the dims based on span
    // create a DMBR for the entries
    unsigned char tmp_DMBR[DMBR_SIZE];
    create_DMBR(split_entries, num_of_entries, tmp_DMBR);

    #ifdef LOG_MBR_SPLITTING_LEAF
    logO.log2File("leaf before ");Node::log_DMBR(tmp_DMBR);
    #endif

    // calculate the span of each dimension 
    double span_list[DIM];
    #ifdef LOG_SPLITTING_VECTOR_INDEX
    logO.log2File("-------------------------\n");
    logO.log2File((int)vector_index);logO.log2File("\n");
    #endif

    for(i = 0; i < DIM; i++){
        int tmp_sum = 0;
        for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
        span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    int sorted_dim_list[DIM]; // hold the sorted dimension list
    // sort the dimension in descending order using span
    NDT_stable_sort(DIM, span_list, sorted_dim_list, true);

	//for(i= 0; i < DIM; i++)
	//	sorted_dim_list[i]=i;
	
    int firstIndexWithMoreThan1Letters;

    for(i = 0; i < DIM; i++)
    {
        if(span_list[i]>static_cast<double>(1.0) / alphabet_sizes[sorted_dim_list[i]])//if this dim has more than 1 letter
        {
            firstIndexWithMoreThan1Letters=i;
            break;
        }
    }

    if(firstIndexWithMoreThan1Letters==DIM)
        firstIndexWithMoreThan1Letters=DIM-1;

    //vector<int> shortestDim;
    //shortestDim.push_back(sorted_dim_list[firstIndexWithMoreThan1Letters]);
    //for(i = firstIndexWithMoreThan1Letters+1; i < DIM; i++)
    //{
    //    if(span_list[i]==span_list[firstIndexWithMoreThan1Letters])
    //        shortestDim.push_back(sorted_dim_list[i]);
    //}

    int best_util_dim;
    int most_entries=0;
	int most_letters=0;
	float best_span;
    //for(i=0;i<shortestDim.size();i++)
    //{
   for(i = firstIndexWithMoreThan1Letters; i < DIM; i++){
  	  int cur_dim=sorted_dim_list[i];
          int num_left_entries;
  		int num_left_letters;
  		if(enforce_minUtil_for_exhaustSplit==0)
  			num_left_entries=sort_entries_by_size_RANDOM_UTIL(cur_dim, split_entries, num_of_entries, sorted_entry_list);
  		else
  			num_left_entries=sort_entries_by_size_RANDOM_UTIL_v3(cur_dim, split_entries, num_of_entries, sorted_entry_list,leaf_min_util,num_left_letters);
  			//num_left_entries=sort_entries_by_size_RANDOM_UTIL_enforceUtil(shortestDim.at(cur_dim), split_entries, num_of_entries, sorted_entry_list,leaf_min_util);
  
          //if(num_left_entries>0)
          //{
          //    best_util_dim=shortestDim.at(cur_dim);
          //    most_entries=num_left_entries;
          //    break;
          //}
  
  		bool findBetterOne=false;
  		if((most_letters==0)&&(num_left_letters>0)){
        best_util_dim=cur_dim;
        most_entries=num_left_entries;
  			most_letters=num_left_letters;
  			best_span = span_list[i];
  			findBetterOne=true;
      }
  		else if((most_letters>0)&&(num_left_letters>0)&&(num_left_letters<most_letters)&&(span_list[i]<=best_span))
  		{
              best_util_dim=cur_dim;
              most_entries=num_left_entries;
  			most_letters=num_left_letters;
  			best_span = span_list[i];
  			findBetterOne=true;
  
  		}
  
  		if(findBetterOne)
  		{
  			count1 = most_entries;
  			count2 = num_of_entries - count1;
  			for(j = 0; j < num_of_entries; j++)
  			{
  				if(j < count1)
  					group1[j] = sorted_entry_list[j];
  				else 
  					group2[j - count1] = sorted_entry_list[j];
  			}
  		}
    }
   if (most_letters==0)
    {
        //cout<<"\n resort to old splitting algorithm at leaf level\n"<<endl;
        best_util_dim = split_algorithm_after_knapsack(split_entries, leaf_min_util, alphabet_sizes, group1, count1, group2, count2, DMBR1, DMBR2); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
        return best_util_dim;
    }
    else
    {//splitting dim is best_util_dim
    //# of entries on left side is most_entries

  //      //most_entries=sort_entries_by_size_RANDOM_UTIL(best_util_dim, split_entries, num_of_entries, sorted_entry_list);
		//if(enforce_minUtil_for_exhaustSplit==0)
		//	most_entries=sort_entries_by_size_RANDOM_UTIL(best_util_dim, split_entries, num_of_entries, sorted_entry_list);
		//else
		//	most_entries=sort_entries_by_size_RANDOM_UTIL_v2(best_util_dim, split_entries, num_of_entries, sorted_entry_list,leaf_min_util );
		//	//most_entries=sort_entries_by_size_RANDOM_UTIL_enforceUtil(best_util_dim, split_entries, num_of_entries, sorted_entry_list,leaf_min_util );


        //count1 = most_entries;
        //count2 = num_of_entries - count1;
        //for(j = 0; j < num_of_entries; j++)
        //{
        //    if(j < count1)
        //        group1[j] = sorted_entry_list[j];
        //    else 
        //        group2[j - count1] = sorted_entry_list[j];
        //}

    		for(j=0;j<count1; j++)
    			sorted_entry_list[j]=  group1[j] ;
    		for(j=0;j<count2; j++)
    			sorted_entry_list[count1+j]=  group2[j] ;

        create_DMBR(split_entries, sorted_entry_list, count1, DMBR1); // create DMBR for grp1
        create_DMBR(split_entries, sorted_entry_list + count1, num_of_entries - count1, DMBR2); // create DMBR for grp2
    }


    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("new leaf  ");Node::log_DMBR(DMBR1);
        logO.log2File("old leaf  ");Node::log_DMBR(DMBR2);
    #endif
        
    return best_util_dim;
        
}




//do not sort dimensions
//only try to split on first dim with more than 1 letter
void Leaf_node::split_TO_DEATH_RANDOM_UTIL_v3(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
    int i, j, k;

    int num_of_entries = LEAF_NODE_SIZE + 1;
    int sorted_entry_list[LEAF_NODE_SIZE + 1];

    //int split_start = static_cast<int>(floor(LEAF_NODE_SIZE * leaf_min_util)) - 1; // the index in the entries where possible split point starts
    //if(split_start < 0)
    //    split_start = 0;
    //int split_end = static_cast<int>(ceil(LEAF_NODE_SIZE * (1 - leaf_min_util))); // the index in the entries where possible split point ends
    //if(split_end >= num_of_entries - 1)
    //    split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
    //if(split_end < split_start)
    //    split_end = split_start;
    count1 = 0;
    count2 = 0;

    // sort the dims based on span
    // create a DMBR for the entries
    unsigned char tmp_DMBR[DMBR_SIZE];
    create_DMBR(split_entries, num_of_entries, tmp_DMBR);

    #ifdef LOG_MBR_SPLITTING_LEAF
    logO.log2File("leaf before ");Node::log_DMBR(tmp_DMBR);
    #endif

    // calculate the span of each dimension 
    double span_list[DIM];
    #ifdef LOG_SPLITTING_VECTOR_INDEX
    logO.log2File("-------------------------\n");
    logO.log2File((int)vector_index);logO.log2File("\n");
    #endif

    for(i = 0; i < DIM; i++){
        int tmp_sum = 0;
        for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
        span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    int sorted_dim_list[DIM]; // hold the sorted dimension list
    // sort the dimension in descending order using span
    //NDT_stable_sort(DIM, span_list, sorted_dim_list, true);

	for(i= 0; i < DIM; i++)
		sorted_dim_list[i]=i;
	
    int firstIndexWithMoreThan1Letters;

    for(i = 0; i < DIM; i++)
    {
        if(span_list[i]>static_cast<double>(1.0) / alphabet_sizes[sorted_dim_list[i]])//if this dim has more than 1 letter
        {
            firstIndexWithMoreThan1Letters=i;
            break;
        }
    }

    if(firstIndexWithMoreThan1Letters==DIM)
        firstIndexWithMoreThan1Letters=DIM-1;

    vector<int> shortestDim;
    shortestDim.push_back(sorted_dim_list[firstIndexWithMoreThan1Letters]);
    //for(i = firstIndexWithMoreThan1Letters+1; i < DIM; i++)
    //{
    //    if(span_list[i]==span_list[firstIndexWithMoreThan1Letters])
    //        shortestDim.push_back(sorted_dim_list[i]);
    //}

    int best_util_dim;
    int most_entries=0;
	int most_letters;
    for(i=0;i<shortestDim.size();i++)
    {
        int num_left_entries;
		
		if(enforce_minUtil_for_exhaustSplit==0)
			num_left_entries=sort_entries_by_size_RANDOM_UTIL(shortestDim.at(i), split_entries, num_of_entries, sorted_entry_list);
		else
			num_left_entries=sort_entries_by_size_RANDOM_UTIL_v3(shortestDim.at(i), split_entries, num_of_entries, sorted_entry_list,leaf_min_util,most_letters);
			//num_left_entries=sort_entries_by_size_RANDOM_UTIL_enforceUtil(shortestDim.at(i), split_entries, num_of_entries, sorted_entry_list,leaf_min_util);

        if(num_left_entries>0)
        {
            best_util_dim=shortestDim.at(i);
            most_entries=num_left_entries;
            break;
        }
    }

    if (most_entries==0)
    {
        //cout<<"\n resort to old splitting algorithm at leaf level\n"<<endl;
        split_algorithm_after_knapsack(split_entries, leaf_min_util, alphabet_sizes, group1, count1, group2, count2, DMBR1, DMBR2); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
        return;
    }
    else
    {//splitting dim is best_util_dim
    //# of entries on left side is most_entries

        //most_entries=sort_entries_by_size_RANDOM_UTIL(best_util_dim, split_entries, num_of_entries, sorted_entry_list);
		if(enforce_minUtil_for_exhaustSplit==0)
			most_entries=sort_entries_by_size_RANDOM_UTIL(best_util_dim, split_entries, num_of_entries, sorted_entry_list);
		else
			most_entries=sort_entries_by_size_RANDOM_UTIL_v3(best_util_dim, split_entries, num_of_entries, sorted_entry_list,leaf_min_util ,most_letters);
			//most_entries=sort_entries_by_size_RANDOM_UTIL_enforceUtil(best_util_dim, split_entries, num_of_entries, sorted_entry_list,leaf_min_util );


        count1 = most_entries;
        count2 = num_of_entries - count1;
        for(j = 0; j < num_of_entries; j++)
        {
            if(j < count1)
                group1[j] = sorted_entry_list[j];
            else 
                group2[j - count1] = sorted_entry_list[j];
        }
        create_DMBR(split_entries, sorted_entry_list, count1, DMBR1); // create DMBR for grp1
        create_DMBR(split_entries, sorted_entry_list + count1, num_of_entries - count1, DMBR2); // create DMBR for grp2
    }

    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("new leaf  ");Node::log_DMBR(DMBR1);
        logO.log2File("old leaf  ");Node::log_DMBR(DMBR2);
    #endif
}



void Leaf_node::split_TO_DEATH_RANDOM_UTIL(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
    int i, j, k;

    int num_of_entries = LEAF_NODE_SIZE + 1;
    int sorted_entry_list[LEAF_NODE_SIZE + 1];

    int split_start = static_cast<int>(floor(LEAF_NODE_SIZE * leaf_min_util)) - 1; // the index in the entries where possible split point starts
    if(split_start < 0)
        split_start = 0;
    int split_end = static_cast<int>(ceil(LEAF_NODE_SIZE * (1 - leaf_min_util))); // the index in the entries where possible split point ends
    if(split_end >= num_of_entries - 1)
        split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
    if(split_end < split_start)
        split_end = split_start;
    count1 = 0;
    count2 = 0;

    // sort the dims based on span
    // create a DMBR for the entries
    unsigned char tmp_DMBR[DMBR_SIZE];
    create_DMBR(split_entries, num_of_entries, tmp_DMBR);

    #ifdef LOG_MBR_SPLITTING_LEAF
    logO.log2File("leaf before ");Node::log_DMBR(tmp_DMBR);
    #endif

    // calculate the span of each dimension 
    double span_list[DIM];
    #ifdef LOG_SPLITTING_VECTOR_INDEX
    logO.log2File("-------------------------\n");
    logO.log2File((int)vector_index);logO.log2File("\n");
    #endif

    for(i = 0; i < DIM; i++){
        int tmp_sum = 0;
        for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
        span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    int sorted_dim_list[DIM]; // hold the sorted dimension list
    // sort the dimension in descending order using span
    //NDT_sort(DIM, span_list, sorted_dim_list, false);
    NDT_stable_sort(DIM, span_list, sorted_dim_list, true);


    int firstIndexWithMoreThan1Letters;

    for(i = 0; i < DIM; i++)
    {
        if(span_list[i]>static_cast<double>(1.0) / alphabet_sizes[sorted_dim_list[i]])//if this dim has more than 1 letter
        {
            firstIndexWithMoreThan1Letters=i;
            break;
        }
    }

    if(firstIndexWithMoreThan1Letters==DIM)
        firstIndexWithMoreThan1Letters=DIM-1;

    vector<int> shortestDim;
    shortestDim.push_back(sorted_dim_list[firstIndexWithMoreThan1Letters]);
    for(i = firstIndexWithMoreThan1Letters+1; i < DIM; i++)
    {
        if(span_list[i]==span_list[firstIndexWithMoreThan1Letters])
            shortestDim.push_back(sorted_dim_list[i]);
    }

    int best_util_dim;
    int most_entries=-1;
    for(i=0;i<shortestDim.size();i++)
    {
        int num_left_entries=sort_entries_by_size_RANDOM_UTIL(shortestDim.at(i), split_entries, num_of_entries, sorted_entry_list);
        if(num_left_entries>0)
        {
            best_util_dim=shortestDim.at(i);
            most_entries=num_left_entries;
            break;
        }
    }
    if (most_entries==0)
    {
        //cout<<"\n resort to old splitting algorithm at leaf level\n"<<endl;
        split_algorithm(split_entries, leaf_min_util, alphabet_sizes, group1, count1, group2, count2, DMBR1, DMBR2); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
        return;
    }
    else
    {//splitting dim is best_util_dim
    //# of entries on left side is most_entries

        most_entries=sort_entries_by_size_RANDOM_UTIL(best_util_dim, split_entries, num_of_entries, sorted_entry_list);
        count1 = most_entries;
        count2 = num_of_entries - count1;
        for(j = 0; j < num_of_entries; j++)
        {
            if(j < count1)
                group1[j] = sorted_entry_list[j];
            else 
                group2[j - count1] = sorted_entry_list[j];
        }
        create_DMBR(split_entries, sorted_entry_list, count1, DMBR1); // create DMBR for grp1
        create_DMBR(split_entries, sorted_entry_list + count1, num_of_entries - count1, DMBR2); // create DMBR for grp2
    }

    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("new leaf  ");Node::log_DMBR(DMBR1);
        logO.log2File("old leaf  ");Node::log_DMBR(DMBR2);
    #endif
}

void Leaf_node::split_algorithm_TO_DEATH(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
    int i, j, k;
    int num_of_entries = LEAF_NODE_SIZE + 1;
    int split_start = static_cast<int>(floor(LEAF_NODE_SIZE * leaf_min_util)) - 1; // the index in the entries where possible split point starts
    if(split_start < 0)
        split_start = 0;
    int split_end = static_cast<int>(ceil(LEAF_NODE_SIZE * (1 - leaf_min_util))); // the index in the entries where possible split point ends
    if(split_end >= num_of_entries - 1)
        split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
    if(split_end < split_start)
        split_end = split_start;

    count1 = 0;
    count2 = 0;

    // sort the dims based on span
    // create a DMBR for the entries
    unsigned char tmp_DMBR[DMBR_SIZE];
    create_DMBR(split_entries, num_of_entries, tmp_DMBR);
    #ifdef LOG_MBR_SPLITTING_LEAF
    logO.log2File("leaf before ");Node::log_DMBR(tmp_DMBR);
    #endif
    // calculate the span of each dimension 
    double span_list[DIM];
    #ifdef LOG_SPLITTING_VECTOR_INDEX
    logO.log2File("-------------------------\n");
    logO.log2File((int)vector_index);logO.log2File("\n");
    #endif


    for(i = 0; i < DIM; i++)
    {
        int tmp_sum = 0;
        for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
        span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    int sorted_dim_list[DIM]; // hold the sorted dimension list
    // sort the dimension in descending order using span
    //NDT_sort(DIM, span_list, sorted_dim_list, false);
    NDT_stable_sort(DIM, span_list, sorted_dim_list, true);
    // try every dim in order sorted by span
    double best_overlap, cur_overlap; 
    double best_span, cur_span;
    double best_balance, cur_balance; // the differenece of span of two group
    double best_area, cur_area; // the sum of area of the current best split
    int best_dim, cur_dim; // the dimension where to split
    unsigned char cur_DMBR1[DMBR_SIZE], cur_DMBR2[DMBR_SIZE];
    int tmp_sum1, tmp_sum2;

    best_overlap = DOUBLE_INF; // max possible value
    best_span = 0.0;
    best_balance = DOUBLE_INF; // max possible value
    best_area = DOUBLE_INF; // max possible value
    best_dim = 0; 

    int sorted_entry_list[LEAF_NODE_SIZE + 1];
    for(i = 0; i < DIM; i++)
    {
        cur_dim = sorted_dim_list[i];
        cur_span = span_list[i];
        // sort entries on the dimension 
        sort_entries_by_size(cur_dim, split_entries, num_of_entries, sorted_entry_list);
        int best_j;
        for(j = split_start; j <= split_end; j++)
        {
            create_DMBR(split_entries, sorted_entry_list, j + 1, cur_DMBR1); // create DMBR for grp1
            create_DMBR(split_entries, sorted_entry_list + j + 1, num_of_entries - j - 1, cur_DMBR2); // create DMBR for grp2
            cur_overlap = cal_normed_overlap(cur_DMBR1, cur_DMBR2, alphabet_sizes);
            if((best_overlap - cur_overlap) > 0)
            { //cur_overlap < best_overlap
                best_j = j;
                for(k = 0; k < DMBR_SIZE; k++)
                {
                    DMBR1[k] = cur_DMBR1[k];
                    DMBR2[k] = cur_DMBR2[k];
                }
                best_overlap = cur_overlap;
                best_span = cur_span;
                tmp_sum1 = 0;
                tmp_sum2 = 0;
                for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                {
                    tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                    tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                }
                // use normalized best_balance
                best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                if(best_balance < 0)
                    best_balance = - best_balance;
                best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                best_dim = cur_dim;
            }
            else 
                if(false)
                {
                    if((cur_span - best_span) > DOUBLE_THRESHOLD)
                    { // never happen since dimensions are sorted by span
                    best_j = j;
                    for(k = 0; k < DMBR_SIZE; k++)
                    {
                        DMBR1[k] = cur_DMBR1[k];
                        DMBR2[k] = cur_DMBR2[k];
                    }
                    best_span = cur_span;
                    tmp_sum1 = 0;
                    tmp_sum2 = 0;
                    for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                    {
                        tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                        tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                    }
                    // use normalized best_balance
                    best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                    if(best_balance < 0)
                        best_balance = - best_balance;
                    best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                    best_dim = cur_dim;
                    }
                    else 
                        if(((cur_span - best_span) > 0 ? (cur_span - best_span) : (best_span - cur_span)) <= DOUBLE_THRESHOLD)
                        { // cur_span == best_span
                            tmp_sum1 = 0;
                            tmp_sum2 = 0;
                            for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                            {
                                tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                                tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                            }
                            // use normalized balance
                            cur_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                            if(cur_balance < 0)cur_balance = - cur_balance;
                            if((best_balance - cur_balance) > DOUBLE_THRESHOLD)
                            { // cur_balance < best_balance
                                best_j = j;
                                for(k = 0; k < DMBR_SIZE; k++)
                                {
                                    DMBR1[k] = cur_DMBR1[k];
                                    DMBR2[k] = cur_DMBR2[k];
                                }
                                best_balance = cur_balance;
                                best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                                best_dim = cur_dim;
                            }
                            else 
                                if(((cur_balance - best_balance) > 0 ? (cur_balance - best_balance) : (best_balance - cur_balance)) <= DOUBLE_THRESHOLD)
                                { // cur_balance == best_balance
                                    cur_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                                    if(best_area - cur_area > DOUBLE_THRESHOLD)
                                    { // cur_area < best_area
                                        best_j = j;
                                        for(k = 0; k < DMBR_SIZE; k++)
                                        {
                                            DMBR1[k] = cur_DMBR1[k];
                                            DMBR2[k] = cur_DMBR2[k];
                                        }
                                        best_area = cur_area;
                                    best_dim = cur_dim;
                                    }
                                }
                        }
                }
            }
            // remember group1 and group 2 now
            if(best_dim == cur_dim)
            {
                count1 = best_j + 1;
                count2 = num_of_entries - count1;
                for(k = 0; k < num_of_entries; k++)
                {
                    if(k < count1)group1[k] = sorted_entry_list[k];
                    else group2[k - count1] = sorted_entry_list[k];
                }
            }
        }

    #ifdef LOG_SPLIT_DIM_NUMBER
        logO.log2File("leaf dim:");logO.log2File(best_dim);
        logO.log2File("\n");
    #endif
    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("new leaf  ");Node::log_DMBR(DMBR1);
        logO.log2File("old leaf  ");Node::log_DMBR(DMBR2);
    #endif
    #ifdef LOG_ENTRY_ON_SPLITTING_DIM
        logO.log2File("leaf split on dim ");logO.log2File(best_dim);logO.log2File("\n");
        for(int debugp = 0; debugp < num_of_entries;debugp++)
        {//each entry holds one line in log file
            log_OneEntry_OnOneDscDim(&split_entries[sorted_entry_list[debugp]],sorted_entry_list[debugp],best_dim);
            logO.log2File("\n");
        }
    #endif
}

int Leaf_node::split_algorithm_after_knapsack(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
    int i, j, k;
    int num_of_entries = LEAF_NODE_SIZE + 1;
    int split_start = static_cast<int>(floor(LEAF_NODE_SIZE * leaf_min_util)) - 1; // the index in the entries where possible split point starts
    if(split_start < 0)
        split_start = 0;
    int split_end = static_cast<int>(ceil(LEAF_NODE_SIZE * (1 - leaf_min_util))); // the index in the entries where possible split point ends
    if(split_end >= num_of_entries - 1)
        split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
    if(split_end < split_start)
        split_end = split_start;

    count1 = 0;
    count2 = 0;

    // sort the dims based on span
    // create a DMBR for the entries
    unsigned char tmp_DMBR[DMBR_SIZE];
    create_DMBR(split_entries, num_of_entries, tmp_DMBR);
    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("leaf before ");Node::log_DMBR(tmp_DMBR);
    #endif
    // calculate the span of each dimension 
    double span_list[DIM];
    #ifdef LOG_SPLITTING_VECTOR_INDEX
        logO.log2File("-------------------------\n");
        logO.log2File((int)vector_index);logO.log2File("\n");
    #endif
    for(i = 0; i < DIM; i++)
    {
        int tmp_sum = 0;
        for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
        span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    int sorted_dim_list[DIM]; // hold the sorted dimension list
    // sort the dimension in descending order using span
    NDT_sort(DIM, span_list, sorted_dim_list, false);
    // try every dim in order sorted by span
    double best_overlap, cur_overlap; 
    double best_span, cur_span;
    double best_balance, cur_balance; // the differenece of span of two group
    double best_area, cur_area; // the sum of area of the current best split
    int best_dim, cur_dim; // the dimension where to split
    unsigned char cur_DMBR1[DMBR_SIZE], cur_DMBR2[DMBR_SIZE];
    int tmp_sum1, tmp_sum2;

    best_overlap = DOUBLE_INF; // max possible value
    best_span = 0.0;
    best_balance = DOUBLE_INF; // max possible value
    best_area = DOUBLE_INF; // max possible value
    best_dim = 0; 
	int whichHeuristicUsedLast;

    int sorted_entry_list[LEAF_NODE_SIZE + 1];
    for(i = 0; i < DIM; i++)
    {
        cur_dim = sorted_dim_list[i];
        cur_span = span_list[i];
        //if(best_overlap < DOUBLE_THRESHOLD && (best_span - cur_span) > DOUBLE_THRESHOLD) // best_overlap == 0 & cur_span < best_span
        //    break;
        // sort entries on the dimension 
        sort_entries_by_size(cur_dim, split_entries, num_of_entries, sorted_entry_list);
        int best_j;
        for(j = split_start; j <= split_end; j++)
        {
            create_DMBR(split_entries, sorted_entry_list, j + 1, cur_DMBR1); // create DMBR for grp1
            create_DMBR(split_entries, sorted_entry_list + j + 1, num_of_entries - j - 1, cur_DMBR2); // create DMBR for grp2
            cur_overlap = cal_normed_overlap(cur_DMBR1, cur_DMBR2, alphabet_sizes);
            if((best_overlap - cur_overlap) > DOUBLE_THRESHOLD)
            { //cur_overlap < best_overlap
				whichHeuristicUsedLast=1;
                best_j = j;
                for(k = 0; k < DMBR_SIZE; k++)
                {
                    DMBR1[k] = cur_DMBR1[k];
                    DMBR2[k] = cur_DMBR2[k];
                }
                best_overlap = cur_overlap;
                best_span = cur_span;
                tmp_sum1 = 0;
                tmp_sum2 = 0;
                for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                {
                    tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                    tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                }
                // use normalized best_balance
                best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                if(best_balance < 0)
                    best_balance = - best_balance;
                best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                best_dim = cur_dim;
            }
            else if(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD)
                { // cur_overlap == best_overlap
                    //if((cur_span - best_span) > DOUBLE_THRESHOLD)
					if(false)
                    { // never happen since dimensions are sorted by span
                        best_j = j;
                        for(k = 0; k < DMBR_SIZE; k++)
                        {
                            DMBR1[k] = cur_DMBR1[k];
                            DMBR2[k] = cur_DMBR2[k];
                        }
                        best_span = cur_span;
                        tmp_sum1 = 0;
                        tmp_sum2 = 0;
                        for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                        {
                            tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                            tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                        }
                        // use normalized best_balance
                        best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                        if(best_balance < 0)
                            best_balance = - best_balance;
                        best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                        best_dim = cur_dim;
                    }
                    //else if(((cur_span - best_span) > 0 ? (cur_span - best_span) : (best_span - cur_span)) <= DOUBLE_THRESHOLD)
					else if(true)
                        { // cur_span == best_span
                            tmp_sum1 = 0;
                            tmp_sum2 = 0;
                            for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                            {
                                tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                                tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                            }
                        // use normalized balance
                        cur_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                        if(cur_balance < 0)cur_balance = - cur_balance;
                        //if((best_balance - cur_balance) > DOUBLE_THRESHOLD)
						if(false)
                        { // cur_balance < best_balance
                            best_j = j;
                            for(k = 0; k < DMBR_SIZE; k++)
                            {
                                DMBR1[k] = cur_DMBR1[k];
                                DMBR2[k] = cur_DMBR2[k];
                            }
                            best_balance = cur_balance;
                            best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                            best_dim = cur_dim;
                        }
                        //else if(((cur_balance - best_balance) > 0 ? (cur_balance - best_balance) : (best_balance - cur_balance)) <= DOUBLE_THRESHOLD)
						else if(true)
                            { // cur_balance == best_balance
                                cur_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                                if(best_area - cur_area > DOUBLE_THRESHOLD)
                                { // cur_area < best_area
									whichHeuristicUsedLast=4;
                                    best_j = j;
                                    for(k = 0; k < DMBR_SIZE; k++)
                                    {
                                        DMBR1[k] = cur_DMBR1[k];
                                        DMBR2[k] = cur_DMBR2[k];
                                    }
                                    best_area = cur_area;
                                    best_dim = cur_dim;
                                }
                            }
                        }
                }
            }
            // remember group1 and group 2 now
            if(best_dim == cur_dim)
            {
                count1 = best_j + 1;
                count2 = num_of_entries - count1;
                for(k = 0; k < num_of_entries; k++)
                {
                    if(k < count1)group1[k] = sorted_entry_list[k];
                    else group2[k - count1] = sorted_entry_list[k];
                }
            }
        }

    #ifdef LOG_SPLIT_DIM_NUMBER
        logO.log2File("leaf dim:");logO.log2File(best_dim);
        logO.log2File("\n");
    #endif
    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("new leaf  ");Node::log_DMBR(DMBR1);
        logO.log2File("old leaf  ");Node::log_DMBR(DMBR2);
    #endif



		switch (whichHeuristicUsedLast)
		{
		case 1:
			heuristics_overlap_used_leaf++;
			break;
		case 4:
			heuristics_area_used_leaf++;
			break;
		default:
			break;
		}


    return best_dim;

}


void Leaf_node::calEntriesCharsHistogram(int sort_dim, Leaf_entry* sort_entries, int entriesCount, int& alphabetSize, unsigned char* alphabet2, 
                                         int* letterFreq, vector<vector<int> >& letter_entry_set, int curMaxAlphabetSize){
  
  int i;
  unsigned char ch;
  alphabetSize = 0;
  //for(i = 0; i < MAX_ALPHABET_SIZE; i++) letterFreq[i] = 0; // initialize
  for(i = 0; i < curMaxAlphabetSize; i++) letterFreq[i] = 0; // initialize
  
  // cal histogram of the entries
  for(i = 0; i < entriesCount; i++) {
    ch = sort_entries[i].key[sort_dim];
    letter_entry_set[ch].push_back(i); //[ letterFreq[ch] ] = i;
    letterFreq[ch]++;
    if(letterFreq[ch] == 1){
      alphabet2[alphabetSize] = ch;
      alphabetSize++;
    }
  }
  
  
}

void Leaf_node::NDT_sort_subset(int count, bool asc, vector<unsigned int>& splitDims, vector<unsigned char>& splitDimsCount){
  int position, largest, current,tmp_index;
  unsigned char tmp_weight;
  if(asc){
    
    for(position = count - 1; position > 0; position--){
      // pick the index to the largest element
      largest = 0;
      for(current = 1; current <= position; current++)
        if (splitDimsCount[largest] <= splitDimsCount[current])largest = current;
        
        // swap
        tmp_weight = splitDimsCount[largest];
        splitDimsCount[largest] = splitDimsCount[position];
        splitDimsCount[position] = tmp_weight;
        
        tmp_index = splitDims[largest];
        splitDims[largest] = splitDims[position];
        splitDims[position] = tmp_index;
    }
  }
  else{
    
    for(position = 0; position < count - 1; position++){
      // pick the index to the largest element
      largest = count - 1;
      for(current = position; current <= count - 2; current++)
        if(splitDimsCount[largest] <= splitDimsCount[current])largest = current;
        
        // swap
        tmp_weight = splitDimsCount[largest];
        splitDimsCount[largest] = splitDimsCount[position];
        splitDimsCount[position] = tmp_weight;
        
        tmp_index = splitDims[largest];
        splitDims[largest] = splitDims[position];
        splitDims[position] = tmp_index;
    }
    
  } 
  
}    

void Leaf_node::recursiveUtilCal(int& alphabetSize, int minUtilReq, vector<unsigned char>& varAlphabetFreq, int index, bool* minSplitChars,
bool* curSplitChars, int& curMinUtil, int newMinUtil){
  int nextIndex,j;
  
  curSplitChars[index]= true;
  newMinUtil += varAlphabetFreq[index];
  if(newMinUtil >= minUtilReq){
    if(newMinUtil<curMinUtil){
      curMinUtil = newMinUtil;
      for(j=0; j<alphabetSize; j++)  minSplitChars[j] = curSplitChars[j];
    }
    else{
      curSplitChars[index]= false;
      newMinUtil -= varAlphabetFreq[index];
    }
  }
  else{
    for(nextIndex=index+1; nextIndex<alphabetSize;  nextIndex++){
      recursiveUtilCal(alphabetSize, minUtilReq, varAlphabetFreq, nextIndex, minSplitChars, curSplitChars, curMinUtil, newMinUtil);
      curSplitChars[nextIndex]= false;
    }  
  }
  
  
}

int Leaf_node::greedyMinBalanceSplit(int& alphabetSize, unsigned char* alphabet2, vector<unsigned char>& varAlphabetFreq, bool* minSplitChars, double leaf_min_util,
int curMaxAlphabetSize){
  
  int i,j,k,l, minUtilReq = ceil(leaf_min_util*LEAF_NODE_SIZE), curMinUtil = LEAF_NODE_SIZE, newMinUtil,startPoint;
  bool logOut=false, curSplitChars[curMaxAlphabetSize]; 
  
  for(i=0; i<curMaxAlphabetSize; i++) curSplitChars[i] = false; 
  //minSplitChars
  for(i=0; i<alphabetSize; i++){ //for each sorted char, start a recursive traversal
    newMinUtil = 0; 
    recursiveUtilCal(alphabetSize, minUtilReq, varAlphabetFreq, i, minSplitChars, curSplitChars, curMinUtil, newMinUtil);
    for(j=0; j<alphabetSize; j++) curSplitChars[j] = false;
    
  }

  return curMinUtil;
  
}


int Leaf_node::split_algorithm_sprt(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, 
int &count2, unsigned char* DMBR1, unsigned char* DMBR2){
  
  int i, j, k, tmp_sum, curDimCount, curSplitDim, curDimMinUtil, globalMinUtil = LEAF_NODE_SIZE, best_dim,startIndex;
  int num_of_entries = LEAF_NODE_SIZE + 1, cur_dim, numSplitCandidate;
  bool found, logOutput = false;
  
  //vector<unsigned int> splitDims; 
  //vector<unsigned char> splitDimsCount;
  
  unsigned char tmp_DMBR[DMBR_SIZE],splitDimsCount[DIM],cur_span, minSpan;
  int dimCount[DMBR_SIZE], sorted_dim_list[DIM]; // hold the sorted dimension list
  unsigned int splitDims[DIM];
  double span_list[DIM];
  
  create_DMBR(split_entries, num_of_entries, tmp_DMBR);
  
  for(i = 0; i < DIM; i++){
    int tmp_sum = 0;
    for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
      tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
    dimCount[i] = tmp_sum;
    span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
  }
  
  if(logOutput){for(i = 0; i < DIM; i++){ cout<<dimCount[i]<<" ";}cout<<endl;} 
  
  NDT_sort(DIM, span_list, sorted_dim_list, true);
  for(i = 0; i < DIM; i++){
    splitDims[i] = sorted_dim_list[i];
    splitDimsCount[i] = dimCount[ sorted_dim_list[i] ];
  }  
  
  if(logOutput){ cout<<endl; for(i=0; i<DIM; i++){ cout<<sorted_dim_list[i]<<"::"<<static_cast<double>(splitDimsCount[i])<<"\t"; } cout<<endl; }
  
  unsigned char alphabet2[ALPHA]; // letters that appear in all components
  bool minSplitChars[ALPHA] = {false};
  int alphabetSize;
  int letterFreq[ALPHA];
  vector<vector<int> > letter_entry_set(ALPHA, vector<int>());
  vector<unsigned char>globalMinSplitChars;  //min util char list of first group
  
  
  //find the first dimension with span>1
  for(i = 0; i < DIM; i++){
    if(splitDimsCount[i]>1){
      curSplitDim=splitDims[i];
      minSpan = splitDimsCount[i];
      startIndex = i;
      break;    
    }
  }  
  best_dim = -1;
  
  vector<unsigned int> varLengthAlphabet;
  vector<unsigned char> varAlphabetFreq;
  
  for(i = startIndex; i < DIM; i++){
    
    cur_dim = splitDims[i];
    cur_span = splitDimsCount[i];
    
    if(cur_span>minSpan && best_dim>=0)  { break; }
    
    calEntriesCharsHistogram(cur_dim, split_entries, num_of_entries, alphabetSize, alphabet2, letterFreq, letter_entry_set, ALPHA);
    
    if(logOutput){ cout<<"Dim:"<<cur_dim<<endl; for(j=0; j<alphabetSize; j++){ cout<<static_cast<unsigned int>(alphabet2[j])<<"::"<<letterFreq[ alphabet2[j]  ]<<"\t"; } cout<<endl; }
    
    varLengthAlphabet.clear();
    varAlphabetFreq.clear();
    
    for(j=0; j<alphabetSize; j++){
      varLengthAlphabet.push_back(alphabet2[j]);
      varAlphabetFreq.push_back(letterFreq[ alphabet2[j] ]);
    }
    
    NDT_sort_subset(alphabetSize, false, varLengthAlphabet, varAlphabetFreq);
    
    if(logOutput){ cout<<"After FreqSorting. Dim:"<<cur_dim<<endl; 
      for(j=0; j<alphabetSize; j++){ cout<<static_cast<unsigned int>(varLengthAlphabet[j])<<"::"<<static_cast<unsigned int>(varAlphabetFreq[ j ])<<"\t"; } cout<<endl;
    }  
    
    
    curDimMinUtil = greedyMinBalanceSplit(alphabetSize, alphabet2, varAlphabetFreq, minSplitChars, leaf_min_util, ALPHA);
    
    if(logOutput) cout<<"curDimMinUtil:"<<curDimMinUtil<<endl;
    
    if(curDimMinUtil<globalMinUtil){
      best_dim = cur_dim;
      globalMinUtil = curDimMinUtil;
      globalMinSplitChars.clear(); //clear previous min util char list of first group
      for(j=0; j<alphabetSize; j++){
        if(minSplitChars[j]){ globalMinSplitChars.push_back(varLengthAlphabet[j]); }
        //globalMinSplitChars[j] = minSplitChars[j];
        minSplitChars[j] = false;
      }
      if(logOutput) cout<<"cur_best_dim:"<<best_dim<<"\t curDimMinUtil:"<<curDimMinUtil<<endl;
    }
  
  } //dimension loop ends 
  
  if(logOutput) cout<<"global best_dim:"<<best_dim<<"\tglobalMinUtil:"<<globalMinUtil<<endl;
  
  count1 = count2 = 0;
  if(best_dim>=0){ //best dim found
    
    int temp = globalMinSplitChars.size();
    
    //if(logOutput){ query_result_file<<"temp:"<<temp<<endl; for(j=0; j<temp; j++){query_result_file<<static_cast<unsigned int>(globalMinSplitChars[j])<<"\t";}query_result_file<<endl; }
    // create group1 and group2 now
    for(j = 0; j < num_of_entries; j++){
      found = false;
      for(k=0; k<temp; k++){
        
        if(split_entries[j].key[best_dim]==globalMinSplitChars[k]){
          group1[count1++] =j;
          found = true;
          break;
        }
        
      }
      if(!found){
        group2[count2++] =j;
      }
    }
    
    assert(count1 == globalMinUtil); 
    assert((count1 + count2) == LEAF_NODE_SIZE+1); 
    
  }
  else{ //no best dim so far 
    cout<<"No best dim found"<<endl;
    exit(1);
  }
  
  //create_DMBR(split_entries, sorted_entry_list, j + 1, cur_DMBR1); // create DMBR for grp1
  
  create_DMBR(split_entries, group1, count1, DMBR1); // create DMBR for grp1
  create_DMBR(split_entries, group2, count2, DMBR2); // create DMBR for grp2
  if(logOutput){
    cout<<endl;
    printDMBR(DMBR1);
    printDMBR(DMBR2);
  }
  return best_dim;
  
  
}  
  

void Leaf_node::split_algorithm(Leaf_entry* split_entries, double leaf_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
    int i, j, k;
    int num_of_entries = LEAF_NODE_SIZE + 1;
    int split_start = static_cast<int>(floor(LEAF_NODE_SIZE * leaf_min_util)) - 1; // the index in the entries where possible split point starts
    if(split_start < 0)
        split_start = 0;
    int split_end = static_cast<int>(ceil(LEAF_NODE_SIZE * (1 - leaf_min_util))); // the index in the entries where possible split point ends
    if(split_end >= num_of_entries - 1)
        split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
    if(split_end < split_start)
        split_end = split_start;

    count1 = 0;
    count2 = 0;

    // sort the dims based on span
    // create a DMBR for the entries
    unsigned char tmp_DMBR[DMBR_SIZE];
    create_DMBR(split_entries, num_of_entries, tmp_DMBR);
    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("leaf before ");Node::log_DMBR(tmp_DMBR);
    #endif
    // calculate the span of each dimension 
    double span_list[DIM];
    #ifdef LOG_SPLITTING_VECTOR_INDEX
        logO.log2File("-------------------------\n");
        logO.log2File((int)vector_index);logO.log2File("\n");
    #endif
    for(i = 0; i < DIM; i++)
    {
        int tmp_sum = 0;
        for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
        span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    int sorted_dim_list[DIM]; // hold the sorted dimension list
    // sort the dimension in descending order using span
    NDT_sort(DIM, span_list, sorted_dim_list, false);
    // try every dim in order sorted by span
    double best_overlap, cur_overlap; 
    double best_span, cur_span;
    double best_balance, cur_balance; // the differenece of span of two group
    double best_area, cur_area; // the sum of area of the current best split
    int best_dim, cur_dim; // the dimension where to split
    unsigned char cur_DMBR1[DMBR_SIZE], cur_DMBR2[DMBR_SIZE];
    int tmp_sum1, tmp_sum2;

    best_overlap = DOUBLE_INF; // max possible value
    best_span = 0.0;
    best_balance = DOUBLE_INF; // max possible value
    best_area = DOUBLE_INF; // max possible value
    best_dim = 0; 

    int sorted_entry_list[LEAF_NODE_SIZE + 1];
    for(i = 0; i < DIM; i++)
    {
        cur_dim = sorted_dim_list[i];
        cur_span = span_list[i];
        if(best_overlap < DOUBLE_THRESHOLD && (best_span - cur_span) > DOUBLE_THRESHOLD) // best_overlap == 0 & cur_span < best_span
            break;
        // sort entries on the dimension 
        sort_entries_by_size(cur_dim, split_entries, num_of_entries, sorted_entry_list);
        int best_j;
        for(j = split_start; j <= split_end; j++)
        {
            create_DMBR(split_entries, sorted_entry_list, j + 1, cur_DMBR1); // create DMBR for grp1
            create_DMBR(split_entries, sorted_entry_list + j + 1, num_of_entries - j - 1, cur_DMBR2); // create DMBR for grp2
            cur_overlap = cal_normed_overlap(cur_DMBR1, cur_DMBR2, alphabet_sizes);
            if((best_overlap - cur_overlap) > DOUBLE_THRESHOLD)
            { //cur_overlap < best_overlap
                best_j = j;
                for(k = 0; k < DMBR_SIZE; k++)
                {
                    DMBR1[k] = cur_DMBR1[k];
                    DMBR2[k] = cur_DMBR2[k];
                }
                best_overlap = cur_overlap;
                best_span = cur_span;
                tmp_sum1 = 0;
                tmp_sum2 = 0;
                for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                {
                    tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                    tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                }
                // use normalized best_balance
                best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                if(best_balance < 0)
                    best_balance = - best_balance;
                best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                best_dim = cur_dim;
            }
            else
                if(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD)
                { // cur_overlap == best_overlap
                    if((cur_span - best_span) > DOUBLE_THRESHOLD)
                    { // never happen since dimensions are sorted by span
                        best_j = j;
                        for(k = 0; k < DMBR_SIZE; k++)
                        {
                            DMBR1[k] = cur_DMBR1[k];
                            DMBR2[k] = cur_DMBR2[k];
                        }
                        best_span = cur_span;
                        tmp_sum1 = 0;
                        tmp_sum2 = 0;
                        for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                        {
                            tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                            tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                        }
                        // use normalized best_balance
                        best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                        if(best_balance < 0)
                            best_balance = - best_balance;
                        best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                        best_dim = cur_dim;
                    }
                    else
                        if(((cur_span - best_span) > 0 ? (cur_span - best_span) : (best_span - cur_span)) <= DOUBLE_THRESHOLD)
                        { // cur_span == best_span
                            tmp_sum1 = 0;
                            tmp_sum2 = 0;
                            for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++)
                            {
                                tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                                tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
                            }
                        // use normalized balance
                        cur_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
                        if(cur_balance < 0)cur_balance = - cur_balance;
                        if((best_balance - cur_balance) > DOUBLE_THRESHOLD)
                        { // cur_balance < best_balance
                            best_j = j;
                            for(k = 0; k < DMBR_SIZE; k++)
                            {
                                DMBR1[k] = cur_DMBR1[k];
                                DMBR2[k] = cur_DMBR2[k];
                            }
                            best_balance = cur_balance;
                            best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                            best_dim = cur_dim;
                        }
                        else
                            if(((cur_balance - best_balance) > 0 ? (cur_balance - best_balance) : (best_balance - cur_balance)) <= DOUBLE_THRESHOLD)
                            { // cur_balance == best_balance
                                cur_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                                if(best_area - cur_area > DOUBLE_THRESHOLD)
                                { // cur_area < best_area
                                    best_j = j;
                                    for(k = 0; k < DMBR_SIZE; k++)
                                    {
                                        DMBR1[k] = cur_DMBR1[k];
                                        DMBR2[k] = cur_DMBR2[k];
                                    }
                                    best_area = cur_area;
                                    best_dim = cur_dim;
                                }
                            }
                        }
                }
            }
            // remember group1 and group 2 now
            if(best_dim == cur_dim)
            {
                count1 = best_j + 1;
                count2 = num_of_entries - count1;
                for(k = 0; k < num_of_entries; k++)
                {
                    if(k < count1)group1[k] = sorted_entry_list[k];
                    else group2[k - count1] = sorted_entry_list[k];
                }
            }
        }

    #ifdef LOG_SPLIT_DIM_NUMBER
        logO.log2File("leaf dim:");logO.log2File(best_dim);
        logO.log2File("\n");
    #endif
    #ifdef LOG_MBR_SPLITTING_LEAF
        logO.log2File("new leaf  ");Node::log_DMBR(DMBR1);
        logO.log2File("old leaf  ");Node::log_DMBR(DMBR2);
    #endif
}


void Leaf_node::create_DMBR(Leaf_entry* leaf_entries, int num_of_leaf_entries, unsigned char* new_DMBR)
{
    int i, j;
    int byte_no, bit_no;

    // Initialize
    for(i = 0; i < DMBR_SIZE; i++)
        new_DMBR[i] = 0;

    // create the DMBR
    for(i = 0; i < num_of_leaf_entries; i++)
    {
        for(j = 0; j < DIM; j++)
        { // union the new data -- enlarge the DMBR
            byte_no = this->DMBR_byte_lut[j][leaf_entries[i].key[j]];
            bit_no = this->DMBR_bit_lut[j][leaf_entries[i].key[j]];
            new_DMBR[byte_no] |= MASKS[bit_no];
        }
    }
}

void Leaf_node::create_DMBR(Leaf_entry* leaf_entries, int* indices_of_leaf_entries_used, int num_of_leaf_entries_used, unsigned char* new_DMBR){
    int i, j;
    int byte_no, bit_no;

    // Initialize
    for(i = 0; i < DMBR_SIZE; i++)new_DMBR[i] = 0;

    // create the DMBR
    int tmp_index;
    for(i = 0; i < num_of_leaf_entries_used; i++)
    {
        tmp_index = indices_of_leaf_entries_used[i];
        for(j = 0; j < DIM; j++)
        { // union the new data -- enlarge the DMBR
            byte_no = this->DMBR_byte_lut[j][leaf_entries[tmp_index].key[j]];
            bit_no = this->DMBR_bit_lut[j][leaf_entries[tmp_index].key[j]];
            new_DMBR[byte_no] |= MASKS[bit_no];
        }
    }
}


//derived from sort_entries_by_size_RANDOM_UTIL_v2
//sort set entries based on freq desending order(to make letters in one node as less as possible)
int Leaf_node::sort_entries_by_size_RANDOM_UTIL_v3(
	int sort_dim, Leaf_entry* sort_entries, 
	int num_of_sort_entries, 
	int* sorted_entry_list,
	double leaf_min_util,
	int & left_letters)
{
    int i, j;
    letter_entry letter_entry_array[MAX_ALPHABET_SIZE];
    for(i = 0; i < MAX_ALPHABET_SIZE; i++)
        letter_entry_array[i].freq = 0; // initialize
    // cal histogram of the entries
    for(i = 0; i < num_of_sort_entries; i++)
    {
        unsigned char ch = sort_entries[i].key[sort_dim];
        //cout<<(char)(ch+'0')<<endl;
        letter_entry_array[ch].entry_set[letter_entry_array[ch].freq] = i;
        letter_entry_array[ch].freq++;
    }

    int sorted_letter_entry_array[MAX_ALPHABET_SIZE];  // an array of indices into the set_entry_array
    double sorted_letter_entry_array_freq[MAX_ALPHABET_SIZE];  // an array of indices into the set_entry_array

	for(i=0;i< MAX_ALPHABET_SIZE; i++)
	{
		sorted_letter_entry_array[i] = i;
		sorted_letter_entry_array_freq[i]= letter_entry_array[i].freq ;
	
	}

   NDT_stable_sort(MAX_ALPHABET_SIZE, sorted_letter_entry_array_freq, sorted_letter_entry_array, false);



    //// populate
    ////for(i = 0; i < MAX_ALPHABET_SIZE; i++)sorted_letter_entry_array[i] = i;
    //int position,tmp_position1=0,tmp_position2=MAX_ALPHABET_SIZE-1;
    //for(position = 0; position <MAX_ALPHABET_SIZE; position++)
    //{
    //    if(letter_entry_array[position].freq >0)
    //    {
    //        sorted_letter_entry_array[tmp_position1] = position;
    //        tmp_position1++;
    //    }
    //    else
    //    {
    //        sorted_letter_entry_array[tmp_position2] = position;
    //        tmp_position2--;
    //    }
    //}

    int x=0;
    vector<int> freq;
    for(i=0;i<MAX_ALPHABET_SIZE;i++)
    {
        int index= sorted_letter_entry_array[i];
        if(letter_entry_array[index].freq>0)
            freq.push_back(letter_entry_array[index].freq);

        for (j=0;j<letter_entry_array[index].freq;j++)
            sorted_entry_list[x++]=letter_entry_array[index].entry_set[j];
    }
    assert(x==num_of_sort_entries);



    if(freq.size()==1)
	{
		left_letters=0;
        return 0;
	}
    else
	{
		int leastNumEntries = int(leaf_min_util*(count+1));
		int actualNumEntries=0;
		int actualLetters = 0;

#ifdef LOG_KANPSACK_VALUE_AND_WEIGHT	
		   logO.log2File("splitting_height ");logO.log2File(splitting_height );logO.log2File(" ");
		   logO.log2File("knapsack weight:value, ");
		   for(int log_i = 0; log_i<freq.size();log_i++)
		   {
		   
			   logO.log2File(" ,");logO.log2File(freq.at(log_i));logO.log2File(" : ");logO.log2File((int)1);
   
		   }
		   
		   logO.log2File("\n");
	

#endif


		for(i=0;i<freq.size();i++)
		{
			actualLetters++;
			actualNumEntries+=freq.at(i);
			if( (actualNumEntries>=leastNumEntries)&&((count+1-actualNumEntries)>=leastNumEntries))
				break;
		}
		if(actualNumEntries==(count+1))
		{
			left_letters=0;
			return 0;
		}
		else
		{
			left_letters=actualLetters;
			return actualNumEntries;
		}
        //return (freq.at(0));
	}
}


//derived from sort_entries_by_size_RANDOM_UTIL
//this version try to enforce a minimum utilization based on # of entries currently in this node
int Leaf_node::sort_entries_by_size_RANDOM_UTIL_v2(
	int sort_dim, Leaf_entry* sort_entries, 
	int num_of_sort_entries, 
	int* sorted_entry_list,
	double leaf_min_util)
{
    int i, j;
    letter_entry letter_entry_array[MAX_ALPHABET_SIZE];
    for(i = 0; i < MAX_ALPHABET_SIZE; i++)
        letter_entry_array[i].freq = 0; // initialize
    // cal histogram of the entries
    for(i = 0; i < num_of_sort_entries; i++)
    {
        unsigned char ch = sort_entries[i].key[sort_dim];
        //cout<<(char)(ch+'0')<<endl;
        letter_entry_array[ch].entry_set[letter_entry_array[ch].freq] = i;
        letter_entry_array[ch].freq++;
    }

    int sorted_letter_entry_array[MAX_ALPHABET_SIZE];  // an array of indices into the set_entry_array

    // populate
    //for(i = 0; i < MAX_ALPHABET_SIZE; i++)sorted_letter_entry_array[i] = i;
    int position,tmp_position1=0,tmp_position2=MAX_ALPHABET_SIZE-1;
    for(position = 0; position <MAX_ALPHABET_SIZE; position++)
    {
        if(letter_entry_array[position].freq >0)
        {
            sorted_letter_entry_array[tmp_position1] = position;
            tmp_position1++;
        }
        else
        {
            sorted_letter_entry_array[tmp_position2] = position;
            tmp_position2--;
        }
    }

    int x=0;
    vector<int> freq;
    for(i=0;i<MAX_ALPHABET_SIZE;i++)
    {
        int index= sorted_letter_entry_array[i];
        if(letter_entry_array[index].freq>0)
            freq.push_back(letter_entry_array[index].freq);

        for (j=0;j<letter_entry_array[index].freq;j++)
            sorted_entry_list[x++]=letter_entry_array[index].entry_set[j];
    }
    assert(x==num_of_sort_entries);
    if(freq.size()==1)
        return 0;
    else
	{
		int leastNumEntries = int(leaf_min_util*(count+1));
		int actualNumEntries=0;
		for(i=0;i<freq.size();i++)
		{
			actualNumEntries+=freq.at(i);
			if( (actualNumEntries>=leastNumEntries)&&((count+1-actualNumEntries)>=leastNumEntries))
				break;
		}
		if(actualNumEntries==(count+1))
			return 0;
		else
			return actualNumEntries;
        //return (freq.at(0));
	}
}

////derived from sort_entries_by_size_RANDOM_UTIL
////this version try to enforce a minimum utilization based on block size
//int Leaf_node::sort_entries_by_size_RANDOM_UTIL_enforceUtil(
//	int sort_dim, Leaf_entry* sort_entries, 
//	int num_of_sort_entries, 
//	int* sorted_entry_list,
//	double leaf_min_util)
//{
//    int i, j;
//    letter_entry letter_entry_array[MAX_ALPHABET_SIZE];
//    for(i = 0; i < MAX_ALPHABET_SIZE; i++)
//        letter_entry_array[i].freq = 0; // initialize
//    // cal histogram of the entries
//    for(i = 0; i < num_of_sort_entries; i++)
//    {
//        unsigned char ch = sort_entries[i].key[sort_dim];
//        //cout<<(char)(ch+'0')<<endl;
//        letter_entry_array[ch].entry_set[letter_entry_array[ch].freq] = i;
//        letter_entry_array[ch].freq++;
//    }
//
//    int sorted_letter_entry_array[MAX_ALPHABET_SIZE];  // an array of indices into the set_entry_array
//
//    // populate
//    //for(i = 0; i < MAX_ALPHABET_SIZE; i++)sorted_letter_entry_array[i] = i;
//    int position,tmp_position1=0,tmp_position2=MAX_ALPHABET_SIZE-1;
//    for(position = 0; position <MAX_ALPHABET_SIZE; position++)
//    {
//        if(letter_entry_array[position].freq >0)
//        {
//            sorted_letter_entry_array[tmp_position1] = position;
//            tmp_position1++;
//        }
//        else
//        {
//            sorted_letter_entry_array[tmp_position2] = position;
//            tmp_position2--;
//        }
//    }
//
//    int x=0;
//	vector<int> freq;
//    vector< vector< int> > entriesGroup;
//    for(i=0;i<MAX_ALPHABET_SIZE;i++)
//    {
//        int index= sorted_letter_entry_array[i];
//        if(letter_entry_array[index].freq>0)
//            freq.push_back(letter_entry_array[index].freq);
//
//		vector<int> tempVec;
//		tempVec.clear();
//        for (j=0;j<letter_entry_array[index].freq;j++)
//		{
//            sorted_entry_list[x++]=letter_entry_array[index].entry_set[j];
//			tempVec.push_back(letter_entry_array[index].entry_set[j]);
//		}
//		entriesGroup.push_back(tempVec);
//    }
//
//	assert(entriesGroup.size()==freq.size());
//    assert(x==num_of_sort_entries);
//    if(freq.size()==1)
//        return 0;
//    else
//	{
//		int leastNumEntries = (int)(leaf_min_util*(count+1));
//		int actualNumEntries=0;
//		for(i=0;i<freq.size();i++)
//		{
//			actualNumEntries+=freq.at(i);
//			if (actualNumEntries>=leastNumEntries)
//				break;
//		}
//		if(actualNumEntries==(count+1))
//			return 0;
//		else
//			return actualNumEntries;
//        //return (freq.at(0));
//	}
//}
//

int Leaf_node::sort_entries_by_size_RANDOM_UTIL(int sort_dim, Leaf_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list){
    int i, j;
    letter_entry letter_entry_array[MAX_ALPHABET_SIZE];
    for(i = 0; i < MAX_ALPHABET_SIZE; i++)
        letter_entry_array[i].freq = 0; // initialize
    // cal histogram of the entries
    for(i = 0; i < num_of_sort_entries; i++)
    {
        unsigned char ch = sort_entries[i].key[sort_dim];
        //cout<<(char)(ch+'0')<<endl;
        letter_entry_array[ch].entry_set[letter_entry_array[ch].freq] = i;
        letter_entry_array[ch].freq++;
    }

    int sorted_letter_entry_array[MAX_ALPHABET_SIZE];  // an array of indices into the set_entry_array

    // populate
    //for(i = 0; i < MAX_ALPHABET_SIZE; i++)sorted_letter_entry_array[i] = i;
    int position,tmp_position1=0,tmp_position2=MAX_ALPHABET_SIZE-1;
    for(position = 0; position <MAX_ALPHABET_SIZE; position++)
    {
        if(letter_entry_array[position].freq >0)
        {
            sorted_letter_entry_array[tmp_position1] = position;
            tmp_position1++;
        }
        else
        {
            sorted_letter_entry_array[tmp_position2] = position;
            tmp_position2--;
        }
    }

    int x=0;
    vector<int> freq;
    for(i=0;i<MAX_ALPHABET_SIZE;i++)
    {
        int index= sorted_letter_entry_array[i];
        if(letter_entry_array[index].freq>0)
            freq.push_back(letter_entry_array[index].freq);

        for (j=0;j<letter_entry_array[index].freq;j++)
            sorted_entry_list[x++]=letter_entry_array[index].entry_set[j];
    }
    assert(x==num_of_sort_entries);
    if(freq.size()==1)
        return 0;
    else
        return (freq.at(0));
}

void Leaf_node::sort_entries_by_size(int sort_dim, Leaf_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list)
{
    int i, j;
    unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components
    int size_of_alphabet2 = 0;
    letter_entry letter_entry_array[MAX_ALPHABET_SIZE];
    for(i = 0; i < MAX_ALPHABET_SIZE; i++)
        letter_entry_array[i].freq = 0; // initialize
    // cal histogram of the entries
    for(i = 0; i < num_of_sort_entries; i++)
    {
        unsigned char ch = sort_entries[i].key[sort_dim];
        letter_entry_array[ch].entry_set[letter_entry_array[ch].freq] = i;
        letter_entry_array[ch].freq++;
        if(letter_entry_array[ch].freq == 1)
        {
            alphabet2[size_of_alphabet2] = ch;
            size_of_alphabet2++;
        }
    }

    // sort letters in alphabet2 by freq
    int unsorted_letter_count = size_of_alphabet2;
    int weight1 = 0, weight2 = 0;
    int left = 0, right = size_of_alphabet2 - 1; // indices
    int largest, current;
    unsigned char tmp_ch;
    while(unsorted_letter_count > 0)
    {
        if(weight1 <= weight2)
        { // next letter should be put to the left side
            largest = left;
            for(current = left + 1; current <= right; current++) // pick the letter with the largest freq
                if(letter_entry_array[alphabet2[current]].freq > letter_entry_array[alphabet2[largest]].freq)
                    largest = current;
            // swap
            tmp_ch = alphabet2[largest];
            alphabet2[largest] = alphabet2[left];
            alphabet2[left] = tmp_ch;

            weight1 += letter_entry_array[alphabet2[left]].freq;
            left++;
        }
        else
        {
            largest = right;
            for(current = left; current <= right - 1; current++) // pick the letter with the largest freq
                if(letter_entry_array[alphabet2[current]].freq > letter_entry_array[alphabet2[largest]].freq)
                   largest = current;
            // swap
            tmp_ch = alphabet2[largest];
            alphabet2[largest] = alphabet2[right];
            alphabet2[right] = tmp_ch;
            weight2 += letter_entry_array[alphabet2[right]].freq;
            right--;
        }
        unsorted_letter_count--;
    }

    int tmp_cnt = 0;
    for(i = 0; i < size_of_alphabet2; i++)
        for(j = 0; j < letter_entry_array[alphabet2[i]].freq; j++)
        {
            sorted_entry_list[tmp_cnt] = letter_entry_array[alphabet2[i]].entry_set[j];
            tmp_cnt++;
        }
    assert(tmp_cnt == num_of_sort_entries);
}

void Leaf_node::logDMBR()
{
    unsigned char	tmp_DMBR[DMBR_SIZE];
    create_DMBR(entries, count, tmp_DMBR);
    logO.log2File("leaf edge ");Node::log_DMBR(tmp_DMBR);
}

void Leaf_node::log2file()
{
    unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components
    for(int i=0;i<DIM;i++)
    {
        for(int j=0;j<MAX_ALPHABET_SIZE;j++)
            alphabet2[j]=MAX_ALPHABET_SIZE;
        for(int j=0;j<count;j++)
            alphabet2[entries[j].key[i]]=entries[j].key[i];
        int temp=A[0];

        for(int j=0;j<MAX_ALPHABET_SIZE;j++)
            if(alphabet2[j]!=MAX_ALPHABET_SIZE)
            {
                logO.log2File(alphabet2[j]);
                logO.log2File(" ");
                temp--;
            }

        for(int j=0;j<temp;j++)
            logO.log2File("  ");
        logO.log2File(";\t");
    }
    logO.log2File("\n");
}

void Leaf_node::log_OneEntry_OnOneDscDim(const Leaf_entry* const oneEntry,int entryIndex, int dimNum)
{
    logO.log2File((int)(((*oneEntry).key[dimNum])));logO.log2File(" ");
}


void Leaf_node::print_OneEntry_OnOneDscDim(const Leaf_entry* const oneEntry,int entryIndex, int dimNum)
{
    cout<<(int)(((*oneEntry).key[dimNum]))<<" ";
}



 //check the size of bytes after compression
//int Leaf_node::getCompressedDiskSize(unsigned char* allEntries_DMBR,int numEntries, int *alphabet_sizes)
//{
//	assert(DIR_NODE_OVERHEAD==4);
//	int totalBytes=DIR_NODE_OVERHEAD+ numEntries*(DMBR_SIZE + sizeof(unsigned int)+usedBytesForEntryHeader);
//
//	int p =0;
//	for (int e =0;e<numEntries;e++)
//	{
//		for(int i =0;i<DIM;i++)
//		{
//			int tmp_sum = 0;
//			for(int j = 0; j <BYTES_PER_DIM_IN_DMBR; j++)
//			{
//				tmp_sum += bit_1_count_lut[allEntries_DMBR[p]];
//				p++;
//			}
//
//			if(tmp_sum==alphabet_sizes[i])
//				totalBytes -= BYTES_PER_DIM_IN_DMBR;
//		}
//	}
//
//	return totalBytes;
//
//}
//
//
