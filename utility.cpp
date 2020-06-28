#ifndef __UTILITY__
#define __UTILITY__


#include <math.h>
#include <iostream>
#include "utility.h"


using namespace std;

int transDecimalToBinary (int decimal) {
    int i, j, target, answer;

    target = decimal;
    answer = 0;

    for (i = 0; target > 0; i++) {
        answer = answer + (target % 2) * pow (10, i);
        target = target / 2;
        cout << decimal << " " << target % 2 << " " << (target % 2) * pow (10, i) << " " << target << " " << answer << endl;
    }
    cout << endl;
    return answer;
}



#endif