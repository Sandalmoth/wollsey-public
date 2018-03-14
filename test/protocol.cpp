#include <iostream>
#include <cassert>
#include <chrono>

#include "../protocol.h"


using namespace std;


int main() {
  auto t1 = [val = chrono::system_clock::now()] { return chrono::system_clock::now() - val; };
  cout << chrono::duration_cast<chrono::milliseconds>(t1()).count() << " ms" << endl;;
}
