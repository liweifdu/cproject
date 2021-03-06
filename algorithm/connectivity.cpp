#include <iostream>
#include <string>

using namespace std;

static const int N = 1000;

int fastsearch(){
    int i, p, q, id[N];
    for (i = 0; i < N; ++i)
        id[i] = i;
    while (cin >> p  >> q){
        int t = id[p];
        if (t == id[q]){
            cout << "the two numbers are connected" << endl;
            cout << id[p] << id[q] << endl;
            continue;
        }
        for (i = 0; i < N; ++i){
            if (id[i] == t)
                id[i] = id[q];
        }
        cout << " " << p << " " << q << endl;
    }
    return 0;
}

int fastmerge(){
    int i, j, p, q, id[N];
    for (i = 0; i < N; ++i)
        id[i] = i;
    while (cin >> p  >> q){
        for (i = p; i != id[i]; i = id[i]);
        for (j = q; j != id[j]; j = id[j]);
        if (i == j){
            cout << "the two numbers are connected" << endl;
            cout << id[p] << id[q] << endl;
            continue;
        }
        id[i] = j;
        cout << " " << p << " " << q << endl;
    }
    return 0;
}

int fastmergeweight(){
    int i, j, p, q, id[N], sz[N];
    for (i = 0; i < N; i++){
        id[i] = i;
        sz[i] = i;
    }
    while (cin >> p >> q){
        for (i = p; i != id[i]; i = id[i]);
        for (j = q; j != id[j]; j = id[j]);
        if (i == j)
            continue;
        if (sz[i] < sz[j]){
            id[i] = j;
            sz[j] += sz[i];
        }
        else {
            id[j] = i;
            sz[i] += sz[j];
        }
        cout << " " << p << " " << q << endl;
    }
    return 0;
}

int equaldivide(){
    int i, j, p, q, id[N], sz[N];
    for (i = 0; i < N; i++){
        id[i] = i;
        sz[i] = i;
    }
    while (cin >> p >> q){
        for (i = p; i != id[i]; i = id[i])
            id[i] = id[id[i]];
        for (j = q; j != id[j]; j = id[j])
            id[j] = id[id[j]];
        if (i == j)
            continue;
        if (sz[i] < sz[j]){
            id[i] = j;
            sz[j] += sz[i];
        }
        else {
            id[j] = i;
            sz[i] += sz[j];
        }
        cout << " " << p << " " << q << endl;
    }
    return 0;
}

int main(){
    //fastsearch();
    //fastmerge();
    //fastmergeweight();
    equaldivide();
}
