#ifndef CONFIG_H
#define CONFIG_H

#include <stack>

// This file defines all the configuration parameters for the ND Tree program

enum ENVIRONMENT {UNIX, WINDOWS};
ENVIRONMENT RUNNING_ENVIRONMENT = UNIX;

const int PartitionType = 2; //1=space partition, 2=data partition


const int QUERY_RESULTS_BUFFER_SIZE = 100000;

typedef int ND_tree_record_count;
typedef int ND_tree_record;
typedef int SplitDimBitBool;

const int MAX_LINE_IN_SOURCE_FILE = 1000000;

const int DIM = 18;
const int CNTDIM = 0; 

//rest 2 are used as part of sourceData file name: sourceData32+10
const int TOTAL_DSC_VALUE =8; //total DSC values to be read out
const int TOTAL_CNT_VALUE =8;//total CNT values to be read out
const int TOTAL_BOX_QUERY_NUM = 100;
const int TOTAL_RANGE_QUERY_NUM =200;

char prDict[] = {'f','l','i','m','v', 's','p','t','a', 'y','*','h','q', 'n','k','d','e', 'c','w','r','g'};
//a = 0, t = 1, g = 3, c = 2                  
string prToDnaDictionary[][3] = {{"1","1","1 2"}, {"1 2","1","0 1 2 3"}, {"0","1","0 1 2"}, {"0","1","3"}, {"3","1","0 1 2 3"}, 
                                 {"0 1","2 3","0 1 2 3"},{"2","2","0 1 2 3"},{"0","2","0 1 2 3"},{"3","2","0 1 2 3"},
                                 {"1","0","1 2"},{"1","0 3","0 3"},{"2","0","1 2"},{"2","0","0 3"},
                                  {"0","0","1 2"},{"0","0","0 3"},{"3","0","1 2"},{"3","0","0 3"},
                                   {"1","3","1 2"},{"1","3","3"},{"0 2","3","0 1 2 3"},{"3","3","0 1 2 3"} };


const int ALPHA=4;
int A[] = {ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA}; 
const int BYTES_PER_DIM_IN_DMBR = ((ALPHA%8)==0)?(ALPHA/8):(ALPHA/8+1);

const int BOX_SIZE_START_AT = 1;
const int BOX_SIZE_STEP = 1;
const int BOX_SIZE_STOP_BEFORE = A[0];
const int RANGE_SIZE_STEP = 1;
const int RANGE_STOP_BEFORE = (DIM + CNTDIM)<10?(DIM + CNTDIM):10;

const int DISK_BLOCK_SIZE = 4096; 
const int DMBR_SIZE = DIM * BYTES_PER_DIM_IN_DMBR; 

int boxSize;

const int usedBytesForEntryHeader = ((DIM%8)==0)?(DIM/8):(DIM/8+1);;
//const int USE_LINK_FOR_BOX_QUERY=1;//don't change, set to 1 in all conditions
//const int USE_LINK_FOR_RANGE_QUERY=1;//don't change, set to 1 in all conditions
		
const int MAX_DIM_AND_CONTDIM = 32; //jan 09, maximum lines for cnt and dsc in boxqueryall file

enum SPLIT_TYPE { ORIGINAL, TO_DEATH/*do not use*/,TO_DEATH_MAX_UTIL/*do not use*/,TO_DEATH_RANDOM_UTIL,SPACE_PARTITION };

// This options decides whether we use NDTree or BoNDTree
//SPLIT_TYPE nodeSplitType  = ORIGINAL;//TO_DEATH_RANDOM_UTIL; 
SPLIT_TYPE nodeSplitType  = TO_DEATH_RANDOM_UTIL;
//SPLIT_TYPE nodeSplitType  = SPACE_PARTITION;  

enum ENUM_TREE_TYPE {STATIC_TREE, DYNAMIC_TREE};ENUM_TREE_TYPE TREE_TYPE = STATIC_TREE;

bool RESORT_2_OLD_SPLIT = true;//if true, might resort to old split when nodeSplitType is TO_DEATH_RANDOM_UTIL

//#define LOG_MBR_SPLITTING_LEAF
//#define LOG_MBR_SPLITTING_DIR

//#define LOG_EDGE_LENGTH_LEAF
//#define LOG_EDGE_LENGTH_DIR

//#define LOG_OLD_AND_NEW_BLOCK_LEAF

////#define LOG_UTILIZATION

//#define LOG_ENTRY_ON_SPLITTING_DIM 
//#define ShorterThanSplitDim

//#define LOG_KANPSACK_VALUE_AND_WEIGHT

int heuristics_overlap_used_leaf =0;
int heuristics_area_used_leaf =0;
int heuristics_overlap_used_dir =0;
int heuristics_area_used_dir =0;
int BestChild_covered_area =0;
int BestChild_notcovered_overlap_enlarge =0;
int BestChild_notcovered_area =0;
int BestChild_notcovered_area_enlarge =0;

const int enforce_minUtil_for_exhaustSplit =1;//always 1
const int try_all_dim_with_minUtil = 1;

int readid_global;
char typeid_global;

string globalAuxFilename = "../data/aux";
string globalRecordFilename = "../data/record";
string globalBQFilename = "../data/box_query_random";

char record_type[9999999][16]={'\0'};

#endif
