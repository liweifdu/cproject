#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

int search(){
    int input, target;
    vector<int> nums;
    cout << "input the target:" <<endl;
    cin >> target;
    cout << "input the array:" << endl;
    while(cin >> input)
        nums.push_back(input);
    auto first = nums.begin(), last = nums.end();
    auto mid = first + (last -first)/2;
    while (mid != last && *mid != target){
        if (target < *mid)
            last = mid;
        else
            first = mid + 1;
        mid = first + (last - first)/2;
    }
    return 0;
}

int search2(){
    int input, target;
    vector<int> nums;
    cout << "input the target:" << endl;
    cin >> target;
    cout << "input the array:" << endl;
    while(cin >> input)
        nums.push_back(input);

    int first = 0, last = nums.size();
    while(first != last){
        const int mid = first + (last - first)/2;
        if (nums[mid] == target)
            cout << mid << endl;
        if (nums[first] <= nums[mid]){
            if(nums[first] <= target && target < nums[mid])
                last = mid;
            else
                first = mid + 1;
        }else{
            if(nums[mid] < target && target <= nums[last-1])
                first = mid + 1;
            else
                last = mid;
        }
    }
    return 0;
}

int stringinput(){
    string input;
    vector<string> nums;
    while (cin >> input){
        nums.push_back(input);
        if(input == "!")
            break;
    }
    for(decltype(nums.size()) i = 0; i != nums.size(); ++i)
            cout << nums[i] << endl;
    return 0;
}

int intinput(){
    int input;
    vector<int> nums;
    while(cin >> input)
        nums.push_back(input);
    for(decltype(nums.size()) i = 0; i < nums.size(); ++i){
        if(nums[i] <= 100)
            nums[i] = nums[i]/10 + 1;
    }

    for(decltype(nums.size()) i = 0; i < nums.size(); ++i){
        cout << nums[i] << endl;
    }
    return 0;
}

int sinput(){
    string s("some string");
    if(s.begin() != s.end()){
        auto it = s.begin();
        *it = toupper(*it);
    }
    cout << s << endl;
    for(auto it = s.begin(); it != s.end() && !isspace(*it); ++it)
        *it = toupper(*it);
    cout << s << endl;
    return 0;
}

int doubleinput(){
    vector<int> nums{1,2,3,4,5,6,7,8,9};
    for(auto it = nums.begin(); it != nums.end(); ++it)
        *it *= 2;
    for(auto it = nums.begin(); it != nums.end(); ++it)
        cout << *it << endl;

    for(auto &i : nums)
        i *= 2;
    for(auto &i : nums)
        cout << i << endl;
    return 0;
}

int array(){
    unsigned scores[11] = {};
    unsigned grade;
    while (cin >> grade){
        if(grade <= 100)
            ++scores[grade/10];
    }
    for(auto i : scores)
        cout << i << " ";
    cout << endl;
    return 0;
}

int arrayassign(){
    int a[10] = {};
    for(size_t i = 0; i < 10; ++i)
        a[i] = i;
    for(size_t i = 0; i < 10; ++i)
        cout << a[i] << endl;
    return 0;
}

int pointarray(){
    int ia[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    //auto ia1(ia);
    //auto ia2(&ia[0]);

    int *pend = end(ia);
    cout << "print way 3:" << endl;
    for (int *p = ia; p != pend; ++p)
        cout << *p << " ";
    cout << endl;

    constexpr size_t sz = 5;
    int arr[sz] = {0, 1, 2, 3, 4};
    int *ip = arr;
    int *ip2 = ip + 4;
    cout << "print way 4:" << endl;
    cout << *ip2 << endl;

    constexpr size_t sz2 = 6;
    int arr2[sz2];
    int *ipend = end(arr2);
    cout << "print way 5:" << endl;
    for(int *p = arr2; p != ipend; ++p)
        *p = 0;
    for(int i = 0; i < 6; ++i)
        cout << arr2[i] << " ";
    cout << endl;

    cout << "print way 6:" << endl;
    vector<int> nums(begin(arr2), end(arr2));
    for(decltype(nums.size()) i = 0; i != nums.size(); ++i)
        cout << nums[i] << " ";
    cout << endl;


    return 0;
}

int multarray(){
    int ia[3][4] = {
        {0, 1, 2, 3},
        {4, 5, 6, 7},
        {8, 9, 10, 11}
    };

    int (&row)[4] = ia[1];
    row[3] = 12;
    cout << "the quote of multi-array:" << endl;
    for(auto &row : ia){
        for(auto col : row)
            cout << col << " ";
        cout << endl;
    }
    
    size_t cnt = 0;
    for(auto &row : ia){
        for(auto &col : row){
            col = cnt;
            ++cnt;
        }
    }

    cout << "the use of point for multi-array:" << endl;
    for(auto p = ia; p != ia + 3; ++p){
        for (auto q = *p; q != *p + 4; ++q)
            cout << *q << ' ';
        cout << endl;
    }

    cout << "stl begin and end used for multi-array:" << endl;
    for(auto p = begin(ia); p != end(ia); ++p){
        for (auto q = begin(*p); q != end(*p); ++q)
            cout << *q << ' ';
        cout << endl;
    }

    return 0;
}



int main(){
    //search();
    //search2();
    //stringinput();
    //intinput();
    //sinput();
    //doubleinput();
    //array();
    //arrayassign();
    //pointarray();
    multarray();

    return 0;
}
