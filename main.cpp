#include <iostream>
#include <ctime>
#include <random>
#include <vector>
#include<string>
#include<fstream>
#include <iomanip>
#include<algorithm>
#include <typeinfo>
using namespace std;

struct Instance
{
    long long int items_number;
    long long int stocks_number;
    vector<double> stocks_lengths;
    vector<double> items_lengths;
    vector<long long int> items_demands;
    vector<double> item_profits;
    vector<long long int> stocks_supplys;
    vector<long long int> items_lengths_integer_part;
    vector<long long int> items_lengths_fraction_part;
    vector<long long int> items_profits_integer_part;
    vector<long long int> items_profits_fraction_part;
    vector<long long int> stocks_lengths_integer_part;
    vector<long long int> stocks_lengths_fraction_part;
    vector<long long int> items_demands_decimal_digit;
    vector<long long int> items_profits_decimal_digit;
    long long int stock_supply_decimal_digit;
    long long int digit;
};

vector<long long int> GenerateRandomIntegers(long long int lower_bound, long long int upper_bound, long long int number) {
    vector<long long int> random_integers;
    for (long long int i = 0; i < number; i++) {
        random_integers.push_back(rand() % (upper_bound - lower_bound + 1) + lower_bound);
    }

    return random_integers;
}

vector<long long int> GenerateRandomIntegersLimitedRange(long long int R,vector<double> ranges, long long int number) {
    vector<long long int> random_integers;
    for (long long int i = 0; i < number; i++) {
        long long int a = rand() % ((long long int)floor(R + ranges[i] / 10) - (long long int)ceil(R - ranges[i] / 10) + 1) + ceil(R - ranges[i] / 10);
        random_integers.push_back(max(a,(long long int)1));
    }

    return random_integers;
}

vector<double> GenerateRationalNumbersWithLimitedDigits(long long int mini_number,long long int lower_bound, long long int upper_bound, int number,
                                                        int digit, vector<long long int>& integers_part, vector<long long int>& fractions_part) {
    integers_part = GenerateRandomIntegers(lower_bound, upper_bound, number);
    fractions_part = GenerateRandomIntegers(0, pow(10,digit)-1, number);
    vector<double> rational_numbers;
    for (int i = 0; i < integers_part.size(); i++) {
        rational_numbers.push_back(integers_part[i] + fractions_part[i] / pow(10, digit));
    }

    return rational_numbers;
}

vector<double> GenerateRationalNumbersWithLimitedDigitsValueRange(long long int number,long long int R,
                                                                  long long int digit, vector<double> ranges, vector<long long int>& integers_part, vector<long long int>& fractions_part) {
    integers_part = GenerateRandomIntegersLimitedRange(R,ranges, number);
    fractions_part = GenerateRandomIntegers(0, pow(10, digit) - 1, number);
    vector<double> rational_numbers;
    for (long long int i = 0; i < integers_part.size(); i++) {
        rational_numbers.push_back(integers_part[i] + fractions_part[i] / pow(10, digit));
    }

    return rational_numbers;
}



/*The CSP's parameters include
N: the number of items (integer)
M: the number of stock types (integer)
l_j: the length of item j, j=1,...,N (rational number)
L_i: the length of stock type i, i=1,...,M (rational number)
d_j: the demand of item j (integer)
D_i: the supply of stock type i (integer)
*/
Instance GenerateInstance(long long int R,long long int items_number,long long int stocks_number, long long int digit, double delta_1, double delta_2) {
    //Generate the lengths of items
    //the value of the length of the item consists of an integer value and a rational value
    //therefore, generating a random rational number with digits limits requires three steps:
    //1) generate randomly integers 2) Generate random integers in the range 0 to 100 and divide these values by 100
    //3) Sum the resulting integers and rational numbers

    Instance instance;
    instance.digit = digit;
    instance.items_number = items_number;
    instance.stocks_number = stocks_number;
    //generate stock types' lengths in the range [5000,10000]
    instance.stocks_lengths = GenerateRationalNumbersWithLimitedDigits(R,5000, 9999, stocks_number, instance.digit,
                                                                       instance.stocks_lengths_integer_part,instance.stocks_lengths_fraction_part);
    //generate items'lengths in the range [delta_1 L_max, delta_2 L_max]
    auto max_position = max_element(instance.stocks_lengths.begin(), instance.stocks_lengths.end());
    double max_stock_length = *max_position;
    instance.items_lengths = GenerateRationalNumbersWithLimitedDigits(0, floor(delta_1*max_stock_length),
                                                                      floor(delta_2 * max_stock_length), items_number, instance.digit,instance.items_lengths_integer_part,instance.items_lengths_fraction_part);
    //generate items' demands in the range [1,100]
    instance.items_demands = GenerateRandomIntegers(1, 5, items_number);
    instance.item_profits= GenerateRationalNumbersWithLimitedDigits(0,floor(0.5*delta_1 * max_stock_length),
                                                                    floor(0.8*delta_2 * max_stock_length), items_number, instance.digit,instance.items_profits_integer_part,instance.items_profits_fraction_part);
    //generate stocks' supply in the range [3,10]
    instance.stocks_supplys = GenerateRandomIntegers(1, 1, stocks_number);

    return instance;
}

void WriteInstance(string file_path, string file_name, Instance instance) {
    ofstream ofs;
    string path = file_path + file_name+".txt";
    ofs.open(path.c_str(), ios::out);
    if (!ofs.is_open()) {
        cout << "CANNOT OPEN" << endl;
    }

    ofs << "items_number" << endl;
    ofs << instance.items_number << endl;
    ofs << "stocks_number" << endl;
    ofs << instance.stocks_number << endl;
    ofs << "digit" << endl;
    ofs << instance.digit << endl;
    ofs << "items" << endl;
    for (long long int i = 0; i < instance.items_number; i++) {
        ofs << std::setprecision(10) <<instance.items_lengths[i] << "\t"<<instance.items_lengths_integer_part[i]<<"\t"<<instance.items_lengths_fraction_part[i]<<"\t"
            << instance.item_profits[i] << "\t" <<instance.items_profits_integer_part[i]<<"\t"<<instance.items_profits_fraction_part[i]<<"\t"
            << instance.items_demands[i] << endl;
    }
    ofs << "stocks" << endl;
    for (long long int i = 0; i < instance.stocks_number; i++) {
        ofs << std::setprecision(10) << instance.stocks_lengths[i] << "\t" <<instance.stocks_lengths_integer_part[i]<<"\t"<<instance.stocks_lengths_fraction_part[i]<<"\t"
            << instance.stocks_supplys[i] << endl;
    }
    ofs.close();
}

vector<double> GenerateKnapsackCapacity(int digit, vector<long long int> item_weights_integer_part, vector<long long int> item_weights_frac_part,
                              vector<long long int> item_demands, vector<long long int>& integer_part, vector<long long int>& frac_part){
    vector<int> a_list;
    vector<int> b_list;
    for (unsigned int i = 0; i < item_weights_integer_part.size(); i++) {
        a_list.push_back(rand() % ((int)floor(item_demands[i]/3) - 0 + 1) + 0);
    }
    for(unsigned int i = 0; i < 7;i++){
        b_list.push_back(rand() % (a_list.size() - 0 + 1) + 0);
    }

    long long int capacity_integer_part = 0;
    long long int capacity_frac_part = 0;

    for(unsigned int i=0;i<b_list.size();i++){
        int index = b_list[i];
        capacity_integer_part += a_list[index]*item_weights_integer_part[index];
        capacity_frac_part += a_list[index]*item_weights_frac_part[index];
        cout<<index<<"_"<<a_list[index]<<"_"<<item_weights_integer_part[index]<<"_"<<item_weights_frac_part[index]<<endl;
    }
    capacity_integer_part+=floor(capacity_frac_part/pow(10,digit));
    capacity_frac_part-=pow(10,digit)*floor(capacity_frac_part/pow(10,digit));

    integer_part.push_back(capacity_integer_part);
    frac_part.push_back(capacity_frac_part);

    return {capacity_integer_part+capacity_frac_part/pow(10,digit)};
}

Instance GenerateUncorrelatedInstances(long long int items_number, long long int R, long long int digit) {
    Instance instance;
    instance.digit = digit;
    instance.items_number = items_number;
    instance.stocks_number = 1;

    instance.items_lengths = GenerateRationalNumbersWithLimitedDigits(0,ceil(R/2), R, instance.items_number, instance.digit, instance.items_lengths_integer_part, instance.items_lengths_fraction_part);
    //generate items' demands in the range [1,100]
    instance.items_demands = GenerateRandomIntegers(1, 5,instance.items_number);
    instance.item_profits = GenerateRationalNumbersWithLimitedDigits(0,ceil(R / 2),R, instance.items_number, instance.digit, instance.items_profits_integer_part, instance.items_profits_fraction_part);
    //generate stocks' supply in the range [3,10]
    instance.stocks_supplys = GenerateRandomIntegers(1, 1, 1);

    //generate stock types' lengths in the range [5000,10000]
    instance.stocks_lengths = GenerateKnapsackCapacity(digit,instance.items_lengths_integer_part,instance.items_lengths_fraction_part,
                                                       instance.items_demands,instance.stocks_lengths_integer_part,instance.stocks_lengths_fraction_part);

    return instance;
}

vector<double> GenerateWeakCorrelatedProfits(long long int digit,long long int R, vector<long long int> item_profits_integer_part,vector<long long int>& integer_parts,
                                             vector<long long int>& frac_parts){
    for(unsigned int i=0;i<item_profits_integer_part.size();i++){
        long long int upper_bound = item_profits_integer_part[i]+R/10;
        long long int lower_bound = item_profits_integer_part[i]-R/10;
        integer_parts.push_back(rand() % (upper_bound - lower_bound + 1) + lower_bound);
    }
    frac_parts = GenerateRandomIntegers(1,  pow(10,digit)-1,item_profits_integer_part.size());

    vector<double> profit_lists;
    for(unsigned int i=0;i<item_profits_integer_part.size();i++){
        profit_lists.push_back(integer_parts[i]+frac_parts[i]/pow(10,digit));
    }
    return profit_lists;
}

Instance GenerateWeakCorrelatedInstances(long long int items_number, long long int R, long long int digit) {
    //Generate the lengths of items
    //the value of the length of the item consists of an integer value and a rational value
    //therefore, generating a random rational number with digits limits requires three steps:
    //1) generate randomly integers 2) Generate random integers in the range 0 to 100 and divide these values by 100
    //3) Sum the resulting integers and rational numbers

    Instance instance;
    instance.digit = digit;
    instance.items_number = items_number;
    instance.stocks_number = 1;

    instance.items_lengths = GenerateRationalNumbersWithLimitedDigits(0,ceil(R/2), R, instance.items_number, instance.digit, instance.items_lengths_integer_part, instance.items_lengths_fraction_part);
    //generate items' demands in the range [1,100]
    instance.items_demands = GenerateRandomIntegers(1, 5,instance.items_number);
    instance.item_profits = GenerateWeakCorrelatedProfits(digit,R,instance.items_lengths_integer_part,
                                                          instance.items_profits_integer_part,instance.items_profits_fraction_part);
    //generate stocks' supply in the range [3,10]
    instance.stocks_supplys = GenerateRandomIntegers(1, 1, 1);

    instance.stocks_lengths = GenerateKnapsackCapacity(digit,instance.items_lengths_integer_part,instance.items_lengths_fraction_part,
                                                       instance.items_demands,instance.stocks_lengths_integer_part,instance.stocks_lengths_fraction_part);

    return instance;
}


vector<double> GenerateStrongCorrelatedProfits(long long int digit, vector<long long int> item_weights_integer_part,
                                               vector<long long int> item_weights_frac_part,vector<long long int>& integer_parts,
                                             vector<long long int>& frac_parts){
    for(unsigned int i=0;i<item_weights_integer_part.size();i++){
        integer_parts.push_back(item_weights_integer_part[i]+10);
    }
    frac_parts = item_weights_frac_part;

    vector<double> profit_lists;
    for(unsigned int i=0;i<item_weights_integer_part.size();i++){
        profit_lists.push_back(integer_parts[i]+frac_parts[i]/pow(10,digit));
    }
    return profit_lists;
}


Instance GenerateStrongCorrelatedInstances(long long int items_number, long long int R, long long int digit) {
    //Generate the lengths of items
    //the value of the length of the item consists of an integer value and a rational value
    //therefore, generating a random rational number with digits limits requires three steps:
    //1) generate randomly integers 2) Generate random integers in the range 0 to 100 and divide these values by 100
    //3) Sum the resulting integers and rational numbers

    Instance instance;
    instance.digit = digit;
    instance.items_number = items_number;
    instance.stocks_number = 1;

    instance.items_lengths = GenerateRationalNumbersWithLimitedDigits(0,ceil(R/2), R, instance.items_number, instance.digit, instance.items_lengths_integer_part, instance.items_lengths_fraction_part);
    //generate items' demands in the range [1,100]
    instance.items_demands = GenerateRandomIntegers(1, 5,instance.items_number);
    instance.item_profits = GenerateStrongCorrelatedProfits(digit,instance.items_lengths_integer_part,instance.items_lengths_fraction_part,
                                                          instance.items_profits_integer_part,instance.items_profits_fraction_part);
    //generate stocks' supply in the range [3,10]
    instance.stocks_supplys = GenerateRandomIntegers(1, 1, 1);

    instance.stocks_lengths = GenerateKnapsackCapacity(digit,instance.items_lengths_integer_part,instance.items_lengths_fraction_part,
                                                       instance.items_demands,instance.stocks_lengths_integer_part,instance.stocks_lengths_fraction_part);

    return instance;
}

Instance GenerateSubsetInstances(long long int items_number, long long int R, long long int digit) {
    //Generate the lengths of items
    //the value of the length of the item consists of an integer value and a rational value
    //therefore, generating a random rational number with digits limits requires three steps:
    //1) generate randomly integers 2) Generate random integers in the range 0 to 100 and divide these values by 100
    //3) Sum the resulting integers and rational numbers

    Instance instance;
    instance.digit = digit;
    instance.items_number = items_number;
    instance.stocks_number = 1;

    instance.items_lengths = GenerateRationalNumbersWithLimitedDigits(0,ceil(R/2), R, instance.items_number, instance.digit, instance.items_lengths_integer_part, instance.items_lengths_fraction_part);
    //generate items' demands in the range [1,100]
    instance.items_demands = GenerateRandomIntegers(1, 5,instance.items_number);
    instance.item_profits = instance.items_lengths;
    instance.items_profits_integer_part=instance.items_lengths_integer_part;
    instance.items_profits_fraction_part = instance.items_lengths_fraction_part;
    //generate stocks' supply in the range [3,10]
    instance.stocks_supplys = GenerateRandomIntegers(1, 1, 1);

    instance.stocks_lengths = GenerateKnapsackCapacity(digit,instance.items_lengths_integer_part,instance.items_lengths_fraction_part,
                                                       instance.items_demands,instance.stocks_lengths_integer_part,instance.stocks_lengths_fraction_part);

    return instance;
}



void UncorrelatedInstances(int index){
    vector<int> N_list =  {10,50,100,150,200,250,300,350,400};
    vector<int> R_list = {100,1000,10000};
    int digit = 2;
    string file_path = "../dataset1/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<N_list.size();j++){
            string file_name  = "unco_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(N_list[j]);
            Instance instance = GenerateUncorrelatedInstances(N_list[j],R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}

void WeakCorrelatedInstances(int index){
    vector<int> N_list =  {10,50,100,150,200,250,300,350,400};
    vector<int> R_list = {100,1000,10000};
    int digit = 2;
    string file_path = "../dataset1/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<N_list.size();j++){
            string file_name  = "weak_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(N_list[j]);
            Instance instance = GenerateWeakCorrelatedInstances(N_list[j],R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}


void StrongCorrelatedInstances(int index){
    vector<int> N_list =  {10,50,100,150,200,250,300,350,400};
    vector<int> R_list = {100,1000,10000};
    int digit = 2;
    string file_path = "../dataset1/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<N_list.size();j++){
            string file_name  = "strong_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(N_list[j]);
            Instance instance = GenerateStrongCorrelatedInstances(N_list[j],R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}

void SubsetCorrelatedInstances(int index){
    vector<int> N_list =  {10,50,100,150,200,250,300,350,400};
    vector<int> R_list = {100,1000,10000};
    int digit = 2;
    string file_path = "../dataset1/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<N_list.size();j++){
            string file_name  = "subset_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(N_list[j]);
            Instance instance = GenerateSubsetInstances(N_list[j],R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}

void FirsDataset(){
    int instance_number = 5;
    for(int i=0;i<instance_number;i++){
        UncorrelatedInstances(i+1);
        WeakCorrelatedInstances(i+1);
        StrongCorrelatedInstances(i+1);
        SubsetCorrelatedInstances(i+1);
    }

}


void UncorrelatedInstancesSecond(int index){
    int N =  200;
    vector<int> K_list =  {2,4,6};
    vector<int> R_list = {100,1000,10000,100000,1000000};
    int digit = 2;
    string file_path = "../dataset2/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<K_list.size();j++){
            string file_name  = "unco_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(K_list[j]);
            Instance instance = GenerateUncorrelatedInstances(N,R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}

void WeakCorrelatedInstancesSecond(int index){
    int N =  200;
    vector<int> K_list =  {2,4,6};
    vector<int> R_list = {100,1000,10000,100000,1000000};
    int digit = 2;
    string file_path = "../dataset2/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<K_list.size();j++){
            string file_name  = "weak_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(K_list[j]);
            Instance instance = GenerateWeakCorrelatedInstances(N,R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}


void StrongCorrelatedInstancesSecond(int index){
    int N =  200;
    vector<int> K_list =  {2,4,6};
    vector<int> R_list = {100,1000,10000,100000,1000000};
    int digit = 2;
    string file_path = "../dataset2/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<K_list.size();j++){
            string file_name  = "strong_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(K_list[j]);
            Instance instance = GenerateStrongCorrelatedInstances(N,R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}

void SubsetCorrelatedInstancesSecond(int index){
    int N =  200;
    vector<int> K_list =  {2,4,6};
    vector<int> R_list = {100,1000,10000,100000,1000000};
    int digit = 2;
    string file_path = "../dataset2/";

    for(unsigned int i=0;i<R_list.size();i++){
        for(unsigned int j=0;j<K_list.size();j++){
            string file_name  = "subset_"+ to_string(index)+"_"+ to_string(R_list[i])+"_"+ to_string(K_list[j]);
            Instance instance = GenerateSubsetInstances(N,R_list[i],digit);
            WriteInstance(file_path,file_name,instance);
        }
    }
}

void SecondDataset(){
    int instance_number = 5;
    for(int i=0;i<instance_number;i++){
        UncorrelatedInstancesSecond(i+1);
        WeakCorrelatedInstancesSecond(i+1);
        StrongCorrelatedInstancesSecond(i+1);
        SubsetCorrelatedInstancesSecond(i+1);
    }
}


void UncorrelatedInstancesThird(int index){
    vector<int> digit_list =  {2,4,6,8,14,16,18};
    int R=1000;
    int N = 200;
    string file_path = "../dataset3/";

    for(unsigned int i=0;i<digit_list.size();i++){
        string file_name  = "unco_"+ to_string(index)+"_"+ to_string(digit_list[i]);
        Instance instance = GenerateUncorrelatedInstances(N,R,digit_list[i]);
        WriteInstance(file_path,file_name,instance);
    }
}

void WeakCorrelatedInstancesThird(int index){
    vector<int> digit_list =  {2,4,6,8,14,16,18};
    int R=1000;
    int N = 200;
    string file_path = "../dataset3/";

    for(unsigned int i=0;i<digit_list.size();i++){
        string file_name  = "weak_"+ to_string(index)+"_"+ to_string(digit_list[i]);
        Instance instance = GenerateWeakCorrelatedInstances(N,R,digit_list[i]);
        WriteInstance(file_path,file_name,instance);
    }
}


void StrongCorrelatedInstancesThird(int index){
    vector<int> digit_list =  {2,4,6,8,14,16,18};
    int R=1000;
    int N = 200;
    string file_path = "../dataset3/";

    for(unsigned int i=0;i<digit_list.size();i++){
        string file_name  = "strong_"+ to_string(index)+"_"+ to_string(digit_list[i]);
        Instance instance = GenerateStrongCorrelatedInstances(N,R,digit_list[i]);
        WriteInstance(file_path,file_name,instance);
    }
}

void SubsetCorrelatedInstancesThird(int index){
    vector<int> digit_list =  {2,4,6,8,14,16,18};
    int R=1000;
    int N = 200;
    string file_path = "../dataset3/";

    for(unsigned int i=0;i<digit_list.size();i++){
        string file_name  = "subset_"+ to_string(index)+"_"+ to_string(digit_list[i]);
        Instance instance = GenerateSubsetInstances(N,R,digit_list[i]);
        WriteInstance(file_path,file_name,instance);
    }
}

void ThirdDataset(){
    int instance_number = 5;
    for(int i=0;i<instance_number;i++){
        UncorrelatedInstancesThird(i+1);
        WeakCorrelatedInstancesThird(i+1);
        StrongCorrelatedInstancesThird(i+1);
        SubsetCorrelatedInstancesThird(i+1);
    }

}

void UncorrelatedInstancesFourth(int index){
    vector<int> digit_list =  {2,4,6};
    vector<long long int> R_list= {1000,10000,100000,1000000};
    int N = 200;
    string file_path = "../dataset4/";

    for(unsigned int i=0;i<R_list.size();i++) {
        for (unsigned int j = 0; j < digit_list.size(); j++) {
            string file_name = "unco_" + to_string(index) + "_"+ to_string(R_list[i])+"_" + to_string(digit_list[j]);
            Instance instance = GenerateUncorrelatedInstances(N, R_list[i], digit_list[j]);
            WriteInstance(file_path, file_name, instance);
        }
    }
}

void WeakCorrelatedInstancesFourth(int index){
    vector<int> digit_list =  {2,4,6};
    vector<long long int> R_list= {1000,10000,100000,1000000};
    int N = 200;
    string file_path = "../dataset4/";

    for(unsigned int i=0;i<R_list.size();i++) {
        for (unsigned int j = 0; j < digit_list.size(); j++) {
            string file_name = "weak_" + to_string(index) + "_"+ to_string(R_list[i])+"_" + to_string(digit_list[j]);
            Instance instance = GenerateWeakCorrelatedInstances(N, R_list[i], digit_list[j]);
            WriteInstance(file_path, file_name, instance);
        }
    }
}


void StrongCorrelatedInstancesFourth(int index){
    vector<int> digit_list =  {2,4,6};
    vector<long long int> R_list= {1000,10000,100000,1000000};
    int N = 200;
    string file_path = "../dataset4/";

    for(unsigned int i=0;i<R_list.size();i++) {
        for (unsigned int j = 0; j < digit_list.size(); j++) {
            string file_name = "strong_" + to_string(index) + "_"+ to_string(R_list[i])+"_" + to_string(digit_list[j]);
            Instance instance = GenerateStrongCorrelatedInstances(N, R_list[i], digit_list[j]);
            WriteInstance(file_path, file_name, instance);
        }
    }
}

void SubsetCorrelatedInstancesFourth(int index){
    vector<int> digit_list =  {2,4,6};
    vector<long long int> R_list= {1000,10000,100000,1000000};
    int N = 200;
    string file_path = "../dataset4/";

    for(unsigned int i=0;i<R_list.size();i++) {
        for (unsigned int j = 0; j < digit_list.size(); j++) {
            string file_name = "subset_" + to_string(index) + "_"+ to_string(R_list[i])+"_" + to_string(digit_list[j]);
            Instance instance = GenerateSubsetInstances(N, R_list[i], digit_list[j]);
            WriteInstance(file_path, file_name, instance);
        }
    }
}

void FourthDataset(){
    int instance_number = 5;
    for(int i=0;i<instance_number;i++){
        UncorrelatedInstancesFourth(i+1);
        WeakCorrelatedInstancesFourth(i+1);
        StrongCorrelatedInstancesFourth(i+1);
        SubsetCorrelatedInstancesFourth(i+1);
    }

}

int main() {
    FirsDataset();
    SecondDataset();
    ThirdDataset();
    FourthDataset();

    return 0;
}

