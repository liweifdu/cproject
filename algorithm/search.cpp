#include <iostream>

using namespace std;

int seqsearch(int *a, int v, int l){
    for(int i = 0; i < l; ++i){
        if(v == a[i])
            //return immediately can stop the code in advance
            cout << "the index of input is " << i << endl;
            return 0;
    }
    return -1;
}

int dichotomy(int *a, int v, int l){
    int afirst = 0, alast = l;
    while(afirst != alast){
    int half = afirst + (alast - afirst)/2;
        if (a[half] == v){
            //return half;
            cout << "the index of input is " << half << endl;
            return 0;
        }
        if (v < a[half])
            alast = half - 1;
        else
            afirst = half + 1;
    }
    return -1;
}

int main(){
    int a[7] = {0, 5 ,7 ,9, 11, 18, 27};
    int l = 7;
    int v = 18;
    dichotomy(a, v, l);
    //seqsearch(a, v, l);
}
