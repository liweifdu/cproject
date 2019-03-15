#include <iostream>
#include <cstring>

static const int N = 1000;
using namespace std;

int main(){
    int a[N];
    memset(a, 1, sizeof(a));
    for(int i = 2; i < N; ++i)
        if(a[i])
            for(int j = i; j*i < N; ++j)
                a[i*j] = 0;
    for(int i = 2; i < N; i++)
        if (a[i])
            cout << " " << i;
    cout << endl;
}

