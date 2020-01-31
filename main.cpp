#include "ND_tree.h"
#include <vector>
#include "logClass.h"
#include <limits.h>
ofstream OutStream;

// GLOBAL VARIABLES - Added by Alok to make interface of this
// tree consistent with the rest of the trees.
string globalQFilename = "../data/q_seq";
string globalRQFilename = "rangequeryAll.txt";
string globalDataFilename = "../data/data_random";
string globalIndexFilename = "ndTree.dat";

//#define MAX_COPIES_OF_DATA_FILE 50

//double dir_min_util = 0.0000000000003;
const double dir_min_util = ((nodeSplitType==ORIGINAL)||(enforce_minUtil_for_exhaustSplit==1))?0.3:0.0000000000003;
const double leaf_min_util = 0.2;

int tmp_DMBR_byte_lut[DIM][MAX_ALPHABET_SIZE]; // Given a dim and a letter, store the coresponding byte in DMBR
int tmp_DMBR_bit_lut[DIM][MAX_ALPHABET_SIZE]; // Given a dim and a letter, store the coresponding bit in DMBR

int debug_boxQ_leaf_accessed=0;
int debug_boxQ_leaf_hit=0;
int debug_boxQ_leaf_hit_peak=0;
vector<int> debug_boxQ_leaf_hit_for_all;

logClass logO;

int debug_height;

int duplicateDataPoints; //used in batchBuild(...), batchGrow(...)

void LocalDMBRInfoCalculation()
{
    int tmp_byte = 0;
    int i, j;
    for(i = 0; i < DIM; i++)
    {
        //tmp_DMBR_start_byte_lut[i] = tmp_byte;

        for(j = 0; j < A[i]; j++)
        {
            tmp_DMBR_byte_lut[i][j] = tmp_byte + j / BITS_PER_BYTE;
            tmp_DMBR_bit_lut[i][j] = j % BITS_PER_BYTE; ////letters must starts from 0
        }


        /** how many bytes occupied by this DBMR[i][j] **/
        tmp_byte += (A[i] - 1) / BITS_PER_BYTE + 1;

        //tmp_DMBR_end_byte_lut[i] = tmp_byte - 1;
    }



}

Dir_entry makeRandomBoxQueryData(ifstream & query_file)
{
    Dir_entry dirEntry;
    for(int j = 0; j < DMBR_SIZE; j++)
        dirEntry.DMBR[j]=0;
    string line;

    for(int j=0;j<DIM;j++)
    {
        int byte_no,bit_no;
        getline(query_file, line);
        istringstream instr(line);
        int v;
        int box_size;
        instr>>box_size;
        for(int t=0;t<box_size;t++)
        {
            instr>>v;
            byte_no = tmp_DMBR_byte_lut[j][v];
            bit_no = tmp_DMBR_bit_lut[j][v];/** todo: here DMBR_bit_lut[j] should support enough letters **/
            dirEntry.DMBR[byte_no] |= MASKS[bit_no];
        }
    }
   
//Why to consume the rest lines. What are the float values? 
/* 	
    for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
        getline(query_file, line);
    //consumes the rest lines holding float values
    for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
        getline(query_file, line);
*/
    return dirEntry;
}

Dir_entry makeBoxQueryData(ifstream & query_file)
{
    Dir_entry dirEntry;

    for(int j = 0; j < DMBR_SIZE; j++)
        dirEntry.DMBR[j]=0;
    string line;

    for(int j=0;j<DIM;j++)
    {
        int byte_no,bit_no;
        getline(query_file, line);
        istringstream instr(line);
        int v;
        for(int t=0;t<boxSize;t++)
        {
            instr>>v;
            byte_no = tmp_DMBR_byte_lut[j][v];
            bit_no = tmp_DMBR_bit_lut[j][v];/** todo: here DMBR_bit_lut[j] should support enough letters **/
            dirEntry.DMBR[byte_no] |= MASKS[bit_no];
        }
    }
    for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
        getline(query_file, line);
    //consumes the rest lines holding float values
    for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
        getline(query_file, line);
    return dirEntry;
}

void batchGrow_with_duplicate(long skipSize, long sizeGrowTo)
{
    ND_tree ndt;
    Leaf_entry new_data/*, query_data, db_data*/;
    Error_code result;
    ndt.read_existing_tree(globalIndexFilename);
    ifstream data_file;
    data_file.open(globalDataFilename.c_str());
    if (data_file.fail())
    {
        cout<<"cant open file "<<globalDataFilename<<endl;
        exit(1);

    }
    int number_of_io = 0;
    long distinctDataPoints=0;
    string line;
    getline(data_file, line);
    for(long i = 0; i < sizeGrowTo; i++)
    {
#ifdef LOG_VECTOR_INDEX
logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif
#ifdef LOG_SPLITTING_VECTOR_INDEX
        vector_index = i;
#endif
        if(i<skipSize)
        {
            getline(data_file, line);
            continue;
        }
        else
        {
            stringstream instr(line);
           for(int j = 0; j < DIM; j++)
           {
               int n;
               instr >> n;
               new_data.key[j] = n - '0';
           }
           new_data.record_count = 1; 
           result = ndt.insert_use_link(new_data, number_of_io);
           if(result == duplicate_error)
            {
                duplicateDataPoints++;
            }
            else
            {
                distinctDataPoints++;
            }
            getline(data_file, line);
        }
    }
    data_file.close();
    cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
    cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
    cout<<"Total data points read:"<<sizeGrowTo <<endl;
    ndt.print_information( );
}

//this one only reads data file up to the give number of data points 
void batchBuild_with_duplicate(long  size)
{
    cout<<"tse"<<endl;
    logO.clearLogs();
    duplicateDataPoints=0;
    ND_tree ndt;
    Leaf_entry new_data/*, query_data, db_data*/;
    Error_code result;
    ndt.create_empty_tree(A, dir_min_util, leaf_min_util, globalIndexFilename);
    ndt.read_existing_tree(globalIndexFilename);
    long num_of_points=size;
    ifstream data_file;
    data_file.open(globalDataFilename.c_str());
    if (data_file.fail())
    {
        cout<<"cant open file "<<globalDataFilename<<endl;
        exit(1);

    }

    int number_of_io = 0;
    string line;
    getline(data_file, line);
    int distinctDataPoints = 0;
    char n;
    for(long  i = 0; i < num_of_points && !data_file.eof(); i++)
    {
#ifdef LOG_VECTOR_INDEX
        logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif
#ifdef LOG_SPLITTING_VECTOR_INDEX
        vector_index = i;
#endif
        istringstream instr(line);
        for(int j = 0; j < DIM; j++)
        {
            instr >> n;
            new_data.key[j] = n - '0';
        }
        new_data.record_count = 1; 
	//new_data.record = 0;
        result = ndt.insert_use_link(new_data, number_of_io);
        if(result == duplicate_error)
        {
            duplicateDataPoints++;
        }
        else
        {
            distinctDataPoints++;
        }
        getline(data_file, line);
        if(i%100000==0) cout<<i<<endl; 
    }
    data_file.close();
    cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
    cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
    cout<<"Total data points read:"<<size <<endl;
    ndt.print_information( );
}


void clear_record()
{
fstream typeid_file;
const char* typeid_filename = (globalRecordFilename+".typeid").c_str();
typeid_file.open(typeid_filename);
if(typeid_file.fail())
{
    cout<<"can't open file "<< typeid_filename <<endl;
    exit(1);
}
typeid_file.clear();
typeid_file.close();

}

void batchGrow_with_duplicate_record_old(long skipSize, long sizeGrowTo)
{
  ND_tree ndt;
  Leaf_entry new_data/*, query_data, db_data*/;
  Error_code result;
  ndt.read_existing_tree(globalIndexFilename);
  ifstream data_file;
  data_file.open(globalDataFilename.c_str());
  if (data_file.fail()){ cout<<"cant open file "<<globalDataFilename<<endl; exit(1); }
  
  int j, number_of_io = 0;
  long distinctDataPoints=0, duplicateDataPoints=0;
  string line;
  char n;
  
  getline(data_file, line);

  for(long i = 0; i < sizeGrowTo; i++){
    if(i<skipSize){
      getline(data_file, line);
      continue;
    }
    else{
      stringstream instr(line);
      cout<<i<<"||"<<line<<endl;
      cout<<"m0:";
      for(j = 0; j < DIM; j++){
        instr >> n;
        new_data.key[j] = n - '0';
      }
      cout<<"m1:";
      new_data.record_count = 1; 
      cout<<"m2:";
      result = ndt.insert_use_link(new_data, number_of_io);
      if(result == duplicate_error){ duplicateDataPoints++; }
      else{ distinctDataPoints++; }
      getline(data_file, line);
      cout<<endl;
    }  
    if(i%100000==0) cout<<i<<endl;
  }  
      
  //..................................
  
  //---------------------------------
  
  data_file.close();
  cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
  cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
  cout<<"Total data points read:"<<sizeGrowTo <<endl;
  ndt.print_information( );

}

void batchGrow_with_duplicate_record(long skipSize, long sizeGrowTo){
  
  ND_tree ndt;
  Leaf_entry new_data/*, query_data, db_data*/;
  Error_code result;
  ndt.read_existing_tree(globalIndexFilename);
  
  ifstream data_file;
  ifstream aux_file;
  data_file.open(globalDataFilename.c_str());
  aux_file.open(globalAuxFilename.c_str());
  if (data_file.fail()){ cout<<"cant open file "<<globalDataFilename<<endl; exit(1); }
  if(aux_file.fail())  { cout<<"cant open file "<<globalAuxFilename<<endl;  exit(1); }
  
  int number_of_io = 0;
  string line;
  string line_aux;
  getline(data_file, line);
  getline(aux_file, line_aux);
  int distinctDataPoints = 0;
  char n;
  char n_aux;
  
  for(long i = 0; i < sizeGrowTo; i++){
#ifdef LOG_VECTOR_INDEX
    logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif
#ifdef LOG_SPLITTING_VECTOR_INDEX
    vector_index = i;
#endif
    if(i<skipSize){
      getline(data_file, line);
      continue;
    }
    else{
      istringstream instr(line);
      istringstream instr_aux(line_aux);
      for(int j = 0; j < DIM; j++){
        instr >> n;
        new_data.key[j] = n - '0';
      }
      
      instr_aux >> readid_global;
      instr_aux >> typeid_global;
      
      new_data.record_count = 1; 
      
      result = ndt.insert_use_link(new_data, number_of_io);
      
      if(result == duplicate_error) duplicateDataPoints++;
      else distinctDataPoints++;
      getline(data_file, line);
      getline(aux_file, line_aux);
      if(i%100000==0) cout<<i<<endl;
      
    }
  }  
  data_file.close();
  cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
  cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
  cout<<"Total data points read:"<<sizeGrowTo <<endl;
  ndt.print_information( );
  
}
  
//Insert kmers with a link to the record(readid,annotation)
void batchBuild_with_duplicate_record(long  size){


    clear_record();

    logO.clearLogs();
    duplicateDataPoints=0;
    ND_tree ndt;
    Leaf_entry new_data/*, query_data, db_data*/;
    Error_code result;
    ndt.create_empty_tree(A, dir_min_util, leaf_min_util, globalIndexFilename);
    ndt.read_existing_tree(globalIndexFilename);
    long num_of_points=size;
    ifstream data_file;
    ifstream aux_file;
    data_file.open(globalDataFilename.c_str());
    aux_file.open(globalAuxFilename.c_str());
    if (data_file.fail())
    {
        cout<<"cant open file "<<globalDataFilename<<endl;
        exit(1);
    }

    if(aux_file.fail())
    {
	cout<<"cant open file "<<globalAuxFilename<<endl;
	exit(1);
    }
    int number_of_io = 0;
    string line;
    string line_aux;
    getline(data_file, line);
    getline(aux_file, line_aux);
    int distinctDataPoints = 0;
    char n;
    char n_aux;
    for(long  i = 0; i < num_of_points && !data_file.eof(); i++)
    {
      //cout<<"i:"<<i<<endl;
#ifdef LOG_VECTOR_INDEX
        logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif
#ifdef LOG_SPLITTING_VECTOR_INDEX
        vector_index = i;
#endif
        istringstream instr(line);
	      istringstream instr_aux(line_aux);
        for(int j = 0; j < DIM; j++)
        {
            instr >> n;
            new_data.key[j] = n - '0';
        }
	
	instr_aux >> readid_global;
	instr_aux >> typeid_global;

        new_data.record_count = 1; 
        
        if(PartitionType==1){
          result = ndt.insert_use_link_sprt(new_data, number_of_io,i);
        }
        else if(PartitionType==2){
          result = ndt.insert_use_link(new_data, number_of_io);
        }
        
        
        if(result == duplicate_error)
        {
            duplicateDataPoints++;
        }
        else
        {
            distinctDataPoints++;
        }
        getline(data_file, line);
	      getline(aux_file, line_aux);
	      if(i%100000==0) cout<<i<<endl;
	    
    }
    data_file.close();
    cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
    cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
    cout<<"Total data points read:"<<size <<endl;
    ndt.print_information( );
}


void batchRangeQuery(    )
{

    Error_code result;

    string input=globalIndexFilename;
    ND_tree ndt;
    Leaf_entry query_data;

    ndt.read_existing_tree(input);


                Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
                int query_results_size;

    string query_fn = globalRQFilename;

    string str_num_of_points;
    int num_of_points;
    ////cout << endl << "How many query records:" << endl;
    ////getline(cin, str_num_of_points);
    ////num_of_points = string_to_int(str_num_of_points);

    num_of_points=TOTAL_RANGE_QUERY_NUM;

    for(int range=0;range<RANGE_STOP_BEFORE;range+=RANGE_SIZE_STEP)
    {
        //start = clock();

        ifstream query_file;
        query_file.open(query_fn.c_str(),ios::in);
        if (query_file.fail())
        {
            cout<<"cant open file "<<query_fn<<endl;
            exit(1);

        }                

        int number_of_io = 0;
        int total_number_of_io = 0;
        int total_results_size = 0;

        string line;
        int n;
        for(int i = 0; i < num_of_points; i++)
        {
            query_results_size = 0;

            //cout <<endl<<"range query for query point "<< i << endl;
            getline(query_file, line);

            if(query_file.bad()||query_file.fail() )

            {
                cout<<"reading query file error\n";
                exit(0);
            }

            istringstream instr(line);



            for(int j = 0; j < DIM; j++)
            {
                instr >> n;
                query_data.key[j] = n;
            }


            result = ndt.range_query_by_Hamming_dist(query_data, range, query_results, query_results_size, number_of_io);

            //cout << "  " << query_results_size << " " << number_of_io << endl;
            total_number_of_io += number_of_io;
            total_results_size += query_results_size;
        }
        query_file.close();

        //finish = clock();
        cout<<total_number_of_io<<" : "<<num_of_points<<endl;
        //elapsed_time = (double)(finish - start) / CLOCKS_PER_SEC;
        //cout << "Time (in seconds): " << elapsed_time / num_of_points << endl;
        cout << "For range: "<<range<<" range query I/O: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
        cout << "Result size: " << static_cast<double>(total_results_size) / num_of_points << endl;




    }//end of for(int range=1;range<RANGE_STOP_BEFORE;++range)


    if(RUNNING_ENVIRONMENT == WINDOWS)
        system("pause");
}

//output the cancer types in the query result to file
void output_records(Leaf_entry* query_results, int query_results_size)
{
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
    unsigned char output[256]={'\0'};
    output[0]=0;
    int type_num;
    int record_id;

for(int k=0; k<query_results_size; k++)
{
    record_id = query_results[k].record;

    typeid_file.seekg(record_id * sizeof(type_array), ios::beg);
    typeid_file.read((char*)type_array, sizeof(type_array));

    type_num = type_array[0];
    int i;

    for(i=1; i<=type_num; i++)
    {
	int p;
	for(p=1; p<=output[0]; p++)
	{
	    if(type_array[i]==output[p])
		break;
	}
	if(p>output[0])
	{
	    output[0]++;
	    output[output[0]]=type_array[i];
	}
    }
}
 
    query_result_file << (unsigned char)(output[0]+'0')<<" ";
    for(int i=1; i<=output[0]; i++)
	query_result_file << output[i]<<" ";
    query_result_file << endl;

    typeid_file.close();
    query_result_file.close();

}


//output the cancer types in the query result to file
void output_records_array(Leaf_entry* query_results, int query_results_size)
{
    fstream typeid_file, query_result_file;
    const char* query_result_filename = (globalBQFilename+ ".result").c_str();
    query_result_file.open(query_result_filename, fstream::in | fstream::out| fstream::app);

    if(query_result_file.fail())
    {
	cout<<"can't open result file "<<query_result_filename<<endl;
    }

    unsigned char type_array[256]={'\0'};
    unsigned char output[256]={'\0'};
    output[0]=0;
    int type_num;
    int record_id;

for(int k=0; k<query_results_size; k++)
{
    record_id = query_results[k].record;
    type_num = record_type[record_id][0];
    int i;

    for(i=1; i<=type_num; i++)
    {
	int p;
	for(p=1; p<=output[0]; p++)
	{
	    if(record_type[record_id][i]==output[p])
		break;
	}
	if(p>output[0])
	{
	    output[0]++;
	    output[output[0]]=record_type[record_id][i];
	}
    }
}
 
    query_result_file << (unsigned char)(output[0]+'0')<<" ";
    for(int i=1; i<=output[0]; i++)
	query_result_file << output[i]<<" ";
    query_result_file << endl;

    typeid_file.close();
    query_result_file.close();

}

void clear_result()
{
fstream  query_result_file;
const char* query_result_filename = (globalBQFilename+".result").c_str();
query_result_file.open(query_result_filename, ios::out | ios::trunc);
if(query_result_file.fail())
{
    cout<<"can't open file "<<query_result_filename<<endl;
    exit(1);
}
query_result_file.clear();
query_result_file.close();
}

void batchRandomBoxQuery_backup(){
  
  clear_result();
  string input=globalIndexFilename;
  ND_tree ndt;
  ndt.read_existing_tree(input);
  LocalDMBRInfoCalculation();
  // box query loop in the original function
  // begins here
  Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
  int query_results_size;
  Dir_entry boxQueryData;
  
  string query_fn=globalBQFilename;
  int num_of_points=TOTAL_BOX_QUERY_NUM;
  ifstream query_file;
  query_file.open(query_fn.c_str());
  if (query_file.fail())
  {
    cout<<"cant open file "<<query_fn<<endl;
    exit(1);
  }
  
  debug_boxQ_leaf_hit_for_all.clear();
  debug_boxQ_leaf_accessed=0;
  
  int number_of_io,dir_io=0 ;
  int total_number_of_io = 0;
  int total_results_size=0;
  int total_record_num = 0;
  int ii=0;
  
  while(!query_file.eof()){
    //cout<<ii<<endl;
    debug_boxQ_leaf_hit_peak=0;
    boxQueryData = makeRandomBoxQueryData(query_file);
    //if(ii>1){
    ndt.box_query(boxQueryData, query_results, query_results_size, number_of_io, dir_io);
    //calculate total record in the result
    for(int k = 0; k < query_results_size; k++)
      total_record_num += query_results[k].record_count;
    
    if(query_results_size>0){
      // output_records(query_results, query_results_size);
      output_records_array(query_results, query_results_size);
    }
    //}
    total_number_of_io += number_of_io;
    //total_results_size += query_results_size;
    //debug_boxQ_leaf_hit_for_all.push_back(debug_boxQ_leaf_hit_peak);
    ii++;
    
  }
  cout<<total_number_of_io<<endl;
  cout<<"ii::"<<ii<<endl;
  //cout<<total_number_of_io<<endl;
  query_file.close();
  
}

Dir_entry makeRandomBoxQueryFromBqArray(vector<string>& boxQueries,int curIndex){
  Dir_entry dirEntry;
  for(int j = 0; j < DMBR_SIZE; j++)
    dirEntry.DMBR[j]=0;
  string line;
  
  for(int j=0;j<DIM;j++)
  {
    int byte_no,bit_no;
    int box_size;
    line = boxQueries.at(curIndex*DIM+j);
    box_size = (line.length()+1)/2;
    istringstream instr(line);
    int v;
    for(int t=0;t<box_size;t++)
    {
      instr>>v;
      byte_no = tmp_DMBR_byte_lut[j][v];
      bit_no = tmp_DMBR_bit_lut[j][v];/** todo: here DMBR_bit_lut[j] should support enough letters **/
      dirEntry.DMBR[byte_no] |= MASKS[bit_no];
    }
  }
  return dirEntry;
  
}

void createKmersFromQseqs(string qfn, int& kmerCount, int& numOFQuerySeq, vector<string>& qKmers){
  
  int j,k,curSeqLength;
  string line,cur_seq,kmer;
  bool newKmer;
  vector<int> dummy;
  
  //box queries from query sequences
  ifstream query_seq_file;
  
  query_seq_file.open(qfn.c_str());
  if (query_seq_file.fail()){ cout<<"cant open file "<<qfn<<endl; exit(1); }
  
  getline(query_seq_file, line);
  
  while(!query_seq_file.eof() ){
  
    if(line.length()>0){
      if(line.at(0)=='>'){
        numOFQuerySeq++;
        cur_seq = "";
      }		
      getline(query_seq_file, line);
    }
    while(line.length()>0 && line.at(0)!='>'){
      cur_seq += line;
      line = "";
      getline(query_seq_file, line);
    }
    for(j=0; j<cur_seq.length(); j++){
      cur_seq.at(j) = tolower(cur_seq.at(j));
    }	
  
    curSeqLength = cur_seq.length() - DIM/3 + 1;
    
    if(cur_seq.length()>0){
      
      for(j=0; j<curSeqLength; j++){
        newKmer = true;
        kmer = cur_seq.substr(j,DIM/3);
        if(kmer.find('x')!=std::string::npos){
          //x value in k-mer, so not a valid k-mer
        }
        else{
          for(k=0; k<kmerCount; k++){
            if(kmer==qKmers.at(k)){
              newKmer = false;
              break;
            }
          }
          if(newKmer){
            qKmers.push_back(kmer);
            kmerCount++;	
          }
          
        }
        
      }
      
    }
    
  
    
  }  
  query_seq_file.close();
  
}

void printBoxQuery(Dir_entry boxQuery){
  for(int j = 0; j < DIM; j++)   cout<<setw(3)<<static_cast<unsigned int>(boxQuery.DMBR[j]);
  cout<<endl;
}


void batchRandomBoxQuery(){
  
    clear_result();
    string input=globalIndexFilename;
    ND_tree ndt;
    ndt.read_existing_tree(input);
    //ndt.print_information();
    LocalDMBRInfoCalculation();
    // box query loop in the original function
    // begins here
    Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
    int query_results_size;
    Dir_entry boxQueryData;

    int i,j,k,n, aminoAcidIndex, numOFQuerySeq=0, kmerCount=0;
    char aminoAcid;	
    vector<string> boxQueries;
    vector<string> qKmers;
    
    createKmersFromQseqs(globalQFilename, kmerCount, numOFQuerySeq, qKmers);
    
    //for(i=0; i<kmerCount; i++){ cout<<qKmers[i]<<endl; }
    cout<<"kmerCount:"<<kmerCount<<endl;
    
    int kLength =  DIM/3;
    for(j=0; j<kmerCount; j++){
      //pKmersStrings[j]="";
      for(k=0; k<kLength; k++){
        aminoAcid = qKmers.at(j).at(k);
        //cout<<aminoAcid<<endl;
        for(n=0; n<sizeof(prDict); n++){
          if(aminoAcid==prDict[n]){
            aminoAcidIndex = n;
            break;
          }
        }
        //pKmersStrings[j] += prToDnaString[aminoAcidIndex];
        for(n=0; n<3; n++){
          //cout<<prToDnaDictionary[aminoAcidIndex][n]<<endl;
          boxQueries.push_back(prToDnaDictionary[aminoAcidIndex][n]);	
        }
      }
    }
    cout<<"Total # of Lines in BoxQueries :"<<boxQueries.size()<<endl;
    
    //exit(1);
    
    string query_fn=globalBQFilename;
    int num_of_points=0;
    ifstream query_file;
    query_file.open(query_fn.c_str());
    if (query_file.fail()){
        cout<<"cant open file "<<query_fn<<endl;
        exit(1);
    }

    debug_boxQ_leaf_hit_for_all.clear();
    debug_boxQ_leaf_accessed=0;

    int number_of_io, dir_io;
    int total_number_of_io = 0, total_dir_io=0;
    int total_results_size=0;
    int total_record_num = 0;
    int ii=0;
   
    for(i=0; i<kmerCount; i++){
    //while(!query_file.eof()){
      
      //if(num_of_points_updated%1000==0){ cout<<num_of_points_updated<<" ";}
     
       // if(i==9){
          //cout<<ii<<endl;
          debug_boxQ_leaf_hit_peak=0;
          //boxQueryData = makeRandomBoxQueryData(query_file);
          boxQueryData = makeRandomBoxQueryFromBqArray(boxQueries,i);
          
          //cout<<qKmers[i]<<endl;
          //printBoxQuery(boxQueryData);
          ndt.box_query(boxQueryData, query_results, query_results_size, number_of_io, dir_io);
        	//calculate total record in the result
        	for(int k = 0; k < query_results_size; k++)
        		total_record_num += query_results[k].record_count;
        	
        	if(query_results_size>0){
        	   //cout<<"query_results_size:"<<query_results_size<<"\t";
        	   // output_records(query_results, query_results_size);
        	   //output_records_array(query_results, query_results_size);
        	   //cout<<endl;
        	}
        
        total_number_of_io += number_of_io;
        total_dir_io += dir_io;
        total_results_size += query_results_size;
        debug_boxQ_leaf_hit_for_all.push_back(debug_boxQ_leaf_hit_peak);
        ii++;
        num_of_points++;
      //}
    
    }


    query_file.close();
    cout<<"boxSize= RANDOM, AvG boxquery I/O: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
    cout<<" AvG matched data point found="<< static_cast<double>(total_results_size)/num_of_points<< endl; 
    cout << " AvG leaf node accessed: " << static_cast<double>(debug_boxQ_leaf_accessed) / num_of_points << endl;
    cout << " AvG Dir node I/O: " << static_cast<double>(total_dir_io) / num_of_points << endl;
    cout<<"total boxquery I/O="<<static_cast<double>(total_number_of_io)<<endl;
    cout<<"total matched data points="<< static_cast<double>(total_results_size)<< endl; 
    cout<<"total matched records = "<< static_cast<double>(total_record_num)<<endl;
   
    //assert(debug_boxQ_leaf_hit_for_all.size()==num_of_points);
    int debug_tatol=0;
    for(unsigned int i=0;i<debug_boxQ_leaf_hit_for_all.size();i++){
        debug_tatol+=debug_boxQ_leaf_hit_for_all.at(i);
    }
    cout << " AvG leaf node hit peak: " << static_cast<double>(debug_tatol) / num_of_points << endl;
    
    if(RUNNING_ENVIRONMENT == WINDOWS) system("pause");
}

void batchBoxQuery()
{
    string input=globalIndexFilename;
    ND_tree ndt;
    ndt.read_existing_tree(input);
    LocalDMBRInfoCalculation();

	for(boxSize=BOX_SIZE_START_AT;boxSize<BOX_SIZE_STOP_BEFORE;boxSize+=BOX_SIZE_STEP)
    {
        //logO.clearLogs();
        Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
        int query_results_size;


        Dir_entry boxQueryData;

        string query_fn="c:\\temp\\boxqueryAll.txt";
        if(RUNNING_ENVIRONMENT == UNIX)
            query_fn="boxqueryAll.txt";
        else
            query_fn="C:\\temp\\boxqueryAll.txt";
        int num_of_points=TOTAL_BOX_QUERY_NUM;

        ifstream query_file;
        query_file.open(query_fn.c_str());
    if (query_file.fail())
    {
        cout<<"cant open file "<<query_fn<<endl;
        exit(1);

    }

        debug_boxQ_leaf_hit_for_all.clear();
        debug_boxQ_leaf_accessed=0;

        int number_of_io, dir_io=0 ;
        int total_number_of_io = 0;

        int total_results_size=0;
        for(int i = 0; i < num_of_points; i++)
        {
            debug_boxQ_leaf_hit_peak=0;
            boxQueryData = makeBoxQueryData(query_file);
            ndt.box_query(boxQueryData, query_results, query_results_size, number_of_io,dir_io);
            total_number_of_io += number_of_io;
            total_results_size += query_results_size;
            debug_boxQ_leaf_hit_for_all.push_back(debug_boxQ_leaf_hit_peak);

        }


        query_file.close();
        cout<<"boxSize= "<<boxSize<< " AvG boxquery I/O: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
        cout<<" AvG matched data point found="<< static_cast<double>(total_results_size)/num_of_points<< endl; 


        cout << " AvG leaf node accessed: " << static_cast<double>(debug_boxQ_leaf_accessed) / num_of_points << endl;

        assert(debug_boxQ_leaf_hit_for_all.size()==num_of_points);
        int debug_tatol=0;
        for(unsigned int i=0;i<debug_boxQ_leaf_hit_for_all.size();i++)
        {
            
            debug_tatol+=debug_boxQ_leaf_hit_for_all.at(i);
        
        }

        cout << " AvG leaf node hit peak: " << static_cast<double>(debug_tatol) / num_of_points << endl;
    }
    if(RUNNING_ENVIRONMENT == WINDOWS)
        system("pause");
}


vector< vector<int> >  makeBoxQueryData_for_linearScan(ifstream & query_file)
{

    //int maxDimAndContDim = 16;//01/09/2007


    string line;

    vector< vector<int> >  boxQueryData;
    boxQueryData.resize(DIM);

    for(int j=0;j<DIM;j++)
    {

        getline(query_file, line);
        istringstream instr(line);
        boxQueryData[j].clear();

        int v;
        for(int t=0;t<(boxSize*A[j])/10;t++)
        {
            instr>>v;
            boxQueryData[j].push_back(v);

        }
    }


    for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
        getline(query_file, line);


    //consumes the rest lines holding float values
    for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
        getline(query_file, line);

    assert(boxQueryData.size()==DIM);
    assert(boxQueryData[0].size()==boxSize);

    return boxQueryData;
}








//result of this LinearScanBoxQuery should be 
void LinearScanBoxQuery(int dataNUM)
{
//    Leaf_entry query_results [QUERY_RESULTS_BUFFER_SIZE];


//    int query_results_size;


    int db_size=dataNUM;
    string query_fn="c:\\temp\\boxqueryAll.txt";


    ifstream query_file;
    query_file.open(query_fn.c_str());
    if (query_file.fail())
    {
        cout<<"cant open file "<<query_fn<<endl;
        exit(1);

    }




    for(boxSize=1;boxSize<=9;boxSize+=1)
    {


        query_file.seekg(0, ios::beg);


        int num_of_points=200;

        int totalNumOfMatches=0;

        for(int i = 0; i < num_of_points; i++) // for every query point
        {

            vector< vector<int> > matchedPoints;
            matchedPoints.clear();

            vector< vector<int> >  boxQueryData;
            boxQueryData = makeBoxQueryData_for_linearScan(query_file);



            vector<int> db_data;

            string db_fn=globalDataFilename;
            ifstream db_file;
            db_file.open(db_fn.c_str());

    if (db_file.fail())
    {
        cout<<"cant open file "<<db_fn<<endl;
        exit(1);

    }




            for(int j = 0; j < db_size; j++)
            {

                //cout<<"query data"<<i<<"db data"<<j<<endl;
                db_data.clear();
                string line;
                getline(db_file, line);
                istringstream dbstr(line);
                for(int k = 0; k < DIM; k++)
                {
                    int n;
                    dbstr >> n;
                    db_data.push_back(n);
                }


                bool matchOnAllDim = true;        
                for(unsigned int s=0;(s<boxQueryData.size())&&matchOnAllDim ;s++)
                {
                    bool matchOnOneDim=false;
                    for(unsigned int x=0;x<boxQueryData[s].size();x++)    
                        if(boxQueryData[s][x]==db_data[s])
                            matchOnOneDim=true;

                    if(!matchOnOneDim)
                        matchOnAllDim=false;

                }



                if(matchOnAllDim)
                {
                    bool findDuplicate=false;
                    for(unsigned int i=0;i<matchedPoints.size();i++)
                        if(matchedPoints.at(i)==db_data)
                            findDuplicate=true;

                    if(!findDuplicate)
                    {
                        matchedPoints.push_back(db_data);
                        //for(int t1=0;t1<boxQueryData.size();t1++)
                        //{
                        //    for(int t2=0;t2<boxQueryData.at(t1).size();t2++)
                        //        cout<<boxQueryData.at(t1).at(t2)<<" ";
                        //    cout<<endl;
 
                        //}    

                        //for(int t2=0;t2<db_data.size();t2++)
                        //    cout<<db_data.at(t2)<<": ";
                        //cout<<endl;


                    }
                }


            } 
            db_file.close();
            totalNumOfMatches+=(int)matchedPoints.size();

        }// end of         for(int i = 0; i < num_of_points; i++) // for every query point




        cout<<"linear scan boxSize="<<boxSize<<endl;
        cout <<  "matched data # from linear scan:" << static_cast<double>(totalNumOfMatches) / num_of_points << endl;

    }
    query_file.close();



}

#define OPT_TREEFILE "-idxfile"
#define OPT_DSCDIM "-dscdim"
#define OPT_LOADFILE "-load_file"
#define OPT_AUXFILE "-aux_file"
#define OPT_RECORDFILE "-record_file"
#define OPT_QFILE "-qfile"
#define OPT_BQFILE "-bqfile"
#define OPT_GINDEX "-gindex"
#define OPT_RQFILE "-rqfile"
#define OPT_RANGE "-range"
#define OPT_SKIP "-skip"
#define OPT_COUNT "-count"
#define OPT_NEW "-newtree"
#define OPT_HELP "-help"
#define UNDEF_STR "__UNDEF__"
#define UNDEF_LONG -999999
struct options
{
    // Name of the data file.
    string datafile;
    string auxfile;
    string recordfile;
    // Name of the index file
    int genindex;
    string qseqfile;	
    string idxfile;
    string bqfile;// box query file
    string rqfile;// range query file
    float range;
    long skip;
    long count;
    int newtree;
    bool help;
};

void get_options(int argc, char **argv, struct options *opt)
{
    string opt_idxfile(OPT_TREEFILE);
    string opt_load_file(OPT_LOADFILE);
    string opt_aux_file(OPT_AUXFILE);
    string opt_record_file(OPT_RECORDFILE);
    string opt_qseqfile(OPT_QFILE);	
    string opt_bqfile(OPT_BQFILE);
    string opt_rqfile(OPT_RQFILE);
    string opt_range(OPT_RANGE);
    string opt_genIndex(OPT_GINDEX);
    string opt_newtree(OPT_NEW);
    string opt_skip(OPT_SKIP);
    string opt_count(OPT_COUNT);
    string opt_help(OPT_HELP);
    for(long i=0;i<argc;i++)
    {
        if(!opt_idxfile.compare(argv[i]))
        {
            opt->idxfile = argv[i+1];
            i++;
        }
        if(!opt_qseqfile.compare(argv[i])){
          opt->qseqfile = argv[i+1];
          i++;
        }
        if(!opt_bqfile.compare(argv[i])){
            opt->bqfile = argv[i+1];
            i++;
        }
        if(!opt_rqfile.compare(argv[i]))
        {
            opt->rqfile = argv[i+1];
            i++;
        }
        if(!opt_load_file.compare(argv[i]))
        {
            opt->datafile = argv[i+1];
            i++;
        }
	if(!opt_aux_file.compare(argv[i]))
	{
	    opt->auxfile = argv[i+1];
	    i++;
	}
	if(!opt_record_file.compare(argv[i]))
	{
	    opt->recordfile = argv[i+1];
	    i++;
	}
        if(!opt_range.compare(argv[i]))
        {
            stringstream s(argv[i+1]);
            s>>(opt->range);
        }
        if(!opt_count.compare(argv[i]))
        {
            stringstream s(argv[i+1]);
            s>>(opt->count);
        }
        if(!opt_genIndex.compare(argv[i]))
        {
          stringstream s(argv[i+1]);
          s>>(opt->genindex);
        }
        if(!opt_skip.compare(argv[i]))
        {
            stringstream s(argv[i+1]);
            s>>(opt->skip);
        }
        if(!opt_newtree.compare(argv[i]))
        {
          stringstream s(argv[i+1]);
          s>>(opt->newtree);
            //opt->newtree = true;
        }
        if(!opt_help.compare(argv[i]))
        {
            opt->help = true;
        }
    }
}

void display_help()
{
    cout<<"This implementation of BoND Tree was modified by Alok"<<endl;
    cout<<"It only works with discrete dimensions and performs random sized box queries."<<endl;
    cout<<"Please feel free to contact Alok at watvealo@cse.msu.edu for any comments"<<endl;
    cout<<"or suggestions"<<endl;
    cout<<"Currently supported options are "<<endl;
    cout<<"1) "<<OPT_TREEFILE<<" : Specifies the template for index filenames. As there"<<endl;
    cout<<"\tare multiple files created, each file will have this name as the prefix"<<endl;
    cout<<"and be appended with the suitable suffix (0,1,2...)"<<endl;
    cout<<"2) "<<OPT_LOADFILE<<" : Specifies the data file to be loaded in the index"<<endl;
    cout<<"\tDatafile must be in CSV format."<<endl;
    cout<<"3) "<<OPT_BQFILE<<" : Name of the box query file"<<endl;
//  cout<<"\t\tBox query is specified as x1,x2,x3:y1,y2:z1"<<endl;
//  cout<<"\t\t\tHere, x1..x3 are characters range along first dimension, y1-y2 is range along second"<<endl;
//  cout<<"dimension and so on"<<endl
    cout<<"4) "<<OPT_RQFILE<<" : Name of the range query file"<<endl;
    cout<<"5) "<<OPT_RANGE<<" : Radius of the range query"<<endl;
    cout<<"6) "<<OPT_SKIP<<" : number of records in the data file that should be skipped"<<endl;
    cout<<"7) "<<OPT_COUNT<<" : Number of data records to be loaded in the database"<<endl;
    cout<<"8) "<<OPT_NEW<<" : Flag indicating that a new index file should be created"<<endl;
    cout<<"\tANY EXISTING INDEX FILE WILL BE DELETED"<<endl;
}

int main(int argc, char *argv[])
{
    if(nodeSplitType==ORIGINAL)
    {
        assert(dir_min_util>0.1);
        TREE_TYPE = STATIC_TREE;
    }
    options opt;
    //Initialize options to Undef values
    // Name of the data file.
    opt.datafile = "../data/data_random";
    opt.auxfile = "../data/aux";
    opt.recordfile = "../data/record";
     //Name of the index file
    opt.idxfile = "../data/index_random";
    opt.qseqfile = "../data/qseq";
    opt.bqfile = "../data/box_query_random";
    opt.rqfile = UNDEF_STR;
    opt.range = 0.;
    opt.genindex = 1;
    opt.skip = 0;
    opt.count = LONG_MAX;
    opt.newtree = 1;
    opt.help = false;
    get_options(argc, argv, &opt);

    if(opt.auxfile.compare(UNDEF_STR))
    {
	cout<<"global aux file name is: "<<opt.auxfile<<endl;
	globalAuxFilename = opt.auxfile;
    }
    if(opt.recordfile.compare(UNDEF_STR))
   {
	cout<<"global aux file name is: "<<opt.recordfile<<endl;
	globalRecordFilename = opt.recordfile;
    }

    if(opt.idxfile.compare(UNDEF_STR))
    {
    //    logO.log2File(opt.idxfile.insert(0, "Index will be created in the file ").c_str());
        cout<<"Index will be created in "<<opt.idxfile<<endl;
        globalIndexFilename = opt.idxfile;
    }
    
    if(opt.genindex > 0){    
      
        if(opt.datafile.compare(UNDEF_STR)){
          globalDataFilename = opt.datafile;
  	      //globalAuxFilename = opt.auxfile;
          //  logO.log2File(opt.datafile.insert(0, "Source data file : ").c_str());
          cout<<"Source data file :"<<opt.datafile<<endl;
          cout<<"Source aux file :"<<opt.auxfile<<endl;
          cout<<"Number of records to load "<<opt.count<<endl;
        }
        cout<<"opt.newtree:"<<opt.newtree<<endl;
        
        if(opt.newtree>0){
            //    log0.log2File("Creatung a new file");
            cout<<"Creating a new index file"<<endl;
            //original method, without link to data record
      	    //batchBuild_with_duplicate(opt.count);
            batchBuild_with_duplicate_record(opt.count);
        }
        else{
            //    log0.log2File("Modifying existing file");
            cout<<"Modifying an existing file "<<endl;
            //batchGrow_with_duplicate(opt.skip, opt.count);
            batchGrow_with_duplicate_record(opt.skip, opt.count);
        }
    }
    
    if(opt.qseqfile.compare(UNDEF_STR)){
        globalQFilename = opt.qseqfile;	
        globalBQFilename = opt.bqfile;
        //logO.log2File(opt.datafile.insert(0, "Box query file : ").c_str());
        cout<<"Box query file "<<opt.bqfile<<endl;
        batchRandomBoxQuery();
        //batchRandomBoxQuery_backup();
    }
    
    if(opt.rqfile.compare(UNDEF_STR))
    {
        globalRQFilename = opt.rqfile;
        //logO.log2File(opt.datafile.insert(0, "Range query file : ").c_str());
        cout<<"Range query file "<<opt.rqfile<<endl;
        batchRangeQuery();
    }
}
