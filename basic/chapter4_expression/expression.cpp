#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

int odd(){
    int input;
    vector<int> nums;
    cout << "input the array:" << endl;
    while(cin >> input)
        nums.push_back(input);
    for(auto it = nums.begin(); it != nums.end(); ++it)
        *it = (*it % 2 != 0) ? (*it * 2) : (*it);
    for(auto &i : nums)
        cout << i << endl;
    return 0;
}








int main(){
    odd();
    return 0;
}
