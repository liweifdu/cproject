#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

using namespace std;

int main(){
    // 2.1 Primitive Built-in Types
    // 2.1.1 Arithmetic type
    // example code
        /* 
        data type decide the operation meaning
        add operation can be used at integer and class
            i = i + j; 
        */

    // exercises
        /* 2.1
            1) size is different
            2) have sign bit or not
            3) size is different(32bit & 64bit), precision size is different (7 & 10) */
        /* 2.2
            float int int*/
    // !!! Note: select data type according to hardware

    // 2.1.2 Type conversions
    // example code
        /*
        data type will auto changed according to the context
            bool b = 42; // b is true
            int i = b; // i equal 1

            unsigned u = 10;
            int i = -42;
            std::cout << i + i << std::endl; // out -84
            std::cout << u + i << std::endl; // out 4294967264 
        */

    // exercises
        /* 2.3
            unsigned u = 10, u2 = 42;
            std::cout << u2 - u << std::endl; // 32
            std::cout << u - u2 << std::endl; // 4294967264

            int i = 10, i2 = 42;
            std::cout << i2 - i << std::endl; // 32
            std::cout << i - i2 << std::endl; // -32
            std::cout << i - u << std::endl;  // 0
            std::cout << u - i << std::endl;  // 0 */
        /* 2.4
            I am correct. */
    // !!! Note: avoid unpredictable or dependent environment operation
    // !!! Note: avoid use unsigned and signed data type mixed

    // 2.1.3 Literals
    // example code
        /* nothing */
    
    // exercises
        /* 2.5
            a) char, wchar_t, string, wstring_t
            b) int, unsigned int, long, unsigned long, 'd10, 'd12
            c) double, float, long double
            d) int, unsigned int, float, float */
        /* 2.6
            int month = 09 is wrong */
        /* 2.7
            a) don't know
            b) 3.14 float
            c) float
            d) long double */
        /* 2.8
            cout << "2M" << endl;
            cout << "2" << "\t" << "M" << endl;  */

    // 2.2 Variables
    // 2.2.1 Variable Definitions
    // example code
        /* 
        define and initialize
            int sum = 0, value,
                units_sold = 0;
            Sales_item item;
            std::string book("0-201-78345-X");

        list initialize
            {
            int units_sold = 0;
            int units_sold = {0};
            int units_sold{0}; // list initialize
            int units_sold(0);
            } // four stetements are all correct

            long double ld = 3.1415926536;
            int a{ld}, b = {ld}; // wrong: list initialize will corrupt because information absent
            int c(ld), d = ld; // correct: conversion and discard some information
        */
    
    // exercises
        /* 2.9
            a) intput an int value;
            b) wrong, list initialize cannot convert type because information absent
            c) correct;
            d) correct but i = 3; */
        /* 2.10
            global_str 0;
            global_int 0;
            local_int undefined;
            local_str undefined; */
    




    // !!! Note: object is a region of memory that has a type
    // !!! Note: Initialization is not assignment. Initialization happens when a variable is given
    //           a value when it is created. Assignment obliterates an object’s current value
    //           and replaces that value with a new one.
    // !!! Note: Uninitialized Variables Cause Run-Time Problems




    return 0;
}
