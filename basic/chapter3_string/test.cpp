#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

using namespace std;

int main(){
   /* string line;
    while(getline(cin,line))
        cout<< line << endl;*/

    /*string word;
    while(cin >> word)
        cout<< word << endl;*/

    /*string word;
    string output;
    while(cin >> word){
        output = output + word + " ";
        if(cin.get() == '\n')
            break;
    }
    cout << output << endl;*/
   
    /*
    string s("i love you!!!");
    decltype(s.size())
        punct_cnt = 0;
    for (auto c : s)
        if(ispunct(c))
            ++punct_cnt;
    cout << punct_cnt << " punctuation characters in " << s << endl;*/

    /*string s("I love you!!!");
    for (auto &c : s)
        c = toupper(c);
    cout << s << endl;*/

    /*
    string s("some string");
    for(decltype(s.size()) index=0; index != s.size() && !isspace(s[index]); ++index)
        s[index] = toupper(s[index]);
    cout << s << endl;*/

    /*const string hexdigits = "0123456789ABCDEF";
    cout << "Enter a series of numbers between 0 and 15" 
         << endl;
    string results;
    size_t n;
    while (cin >> n)
        if (n < hexdigits.size())
            results += hexdigits[n];
    cout << "Your hex numbers is: " << results << endl;*/

    vector<int> nums;
    int input;
    while(cin >> input)
        nums.push_back(input);
    for(decltype(nums.size()) i = 0; i != nums.size(); ++i)
        cout << nums[i] << endl;

    return 0;
}
