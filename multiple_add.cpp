#include <bits/stdc++.h>
#include <benchmarking.hpp>
using namespace std;

void add() {
  int x = 0;
  for (int i = 0; i < 10000000; i++) x += 10;
}

void multiply() {
  int x = 10 * 10000000;
}

int main() {
  TIMING(multiply);
  TIMING(add);
}