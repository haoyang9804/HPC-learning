#pragma once
#include <chrono>
#include <iostream>

#define MARK_FUNCTION std::cout << "Function name: " << __PRETTY_FUNCTION__ << std::endl;

static std::chrono::high_resolution_clock::time_point startTime, endTime;
static double duration;

template<typename R, typename ...Args>
R TIMING(R(*f)(const Args&...), const Args& ...args) {
  startTime = std::chrono::high_resolution_clock::now();
  R r = f(std::forward<Args>(args)...);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by trivial: " << duration << " microseconds." << std::endl;
  return r;
}

template<typename R, typename ...Args>
R TIMING(R(*f)(Args&...), Args& ...args) {
  startTime = std::chrono::high_resolution_clock::now();
  R r = f(args...);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by trivial: " << duration << " microseconds." << std::endl;
  return r;
}

template<typename R, typename ...Args>
R TIMING(R(*f)(Args...), const Args& ...args) {
  startTime = std::chrono::high_resolution_clock::now();
  R r = f(args...);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by trivial: " << duration << " microseconds." << std::endl;
  return r;
}

// template<>
void TIMING(void(*f)()) {
  startTime = std::chrono::high_resolution_clock::now();
  f();
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by trivial: " << duration << " microseconds." << std::endl;
  return;
}