// An attempt of HPC for 2-dimensional convolution
#include <time.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <random>

std::random_device rd;

struct matrix {
  int h, w;
  double** vals;
  bool init = false;
  matrix(int h, int w, bool init = false): h(h), w(w), init(init) {
    if (init) {
      vals = (double**) malloc(sizeof(double*)*h);
      for (int i = 0; i < h; i++)
        vals[i] = (double*) malloc(sizeof(double)*w);
    }
  }
  bool operator<= (const matrix& other) {
    return h <= other.h && w <= other.w;
  }
  double& at(size_t i, size_t j) {
    assert(init);
    return vals[i][j];
  }
  void set(size_t i, size_t j, double val) {
    assert(init);
    assert(i < h && j < w);
    vals[i][j] = val;
  }
  void random_init() {
    assert(init);
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0, 256);
    for (int i = 0; i < h; i++) {
      for (int j = 0; j < w; j++) {
        vals[i][j] = dist(gen);
      }
    }
  }
  ~matrix() {
    if (init) free(vals);
  }

  void print() {
    for (int i = 0; i < h; i++) {
      for (int j = 0; j < w; j++) {
        std::cout << vals[i][j] << "+";
      }
      std::cout << std::endl;
    }
  }

  bool operator==(matrix &m) {
    if (h != m.h || w != m.w) return false;
    for (int i = 0; i < h; i++) {
      for (int j = 0; j < w; j++) {
        if (std::abs(vals[i][j]-m.at(i, j)) > 1e-5) {
          std::cerr << vals[i][j] << " " << m.at(i, j) << std::endl;
          return false;
        }
      }
    }
    return true;
  }

};

// trivial one, no optimization
matrix trivial(matrix& image, matrix& kernel) {
  assert(kernel <= image);
  matrix res(image.h-kernel.h+1, image.w-kernel.w+1, true);
  for (int i = 0; i + kernel.h <= image.h; i++) {
    for (int j = 0; j + kernel.w <= image.w; j++) {
      double avg = 0;
      for (int ii = i; ii < i + kernel.h; ii++) {
        for (int jj = j; jj < j + kernel.w; jj++) {
          avg += image.at(ii, jj);
        }
      }
      avg /= (kernel.h*kernel.w);
      res.set(i, j, avg);
    }
  }
  return res;
}

// memoization with a matrix
matrix memoize(matrix& image, matrix& kernel) {
  assert(kernel <= image);
  matrix res_1(image.h, image.w-kernel.w+1, true);
  matrix res_2(image.h-kernel.h+1, image.w-kernel.w+1, true);
  double sum = 0;

  for (int i = 0; i < image.h; i++) {
    for (int j = 0; j + kernel.w <= image.w; j++) {
      if (j == 0) {
        sum = 0;
        for (int jj = 0; jj < j + kernel.w; jj++) {
          sum += image.at(i, jj);
        }
      }
      else {
        sum = sum - image.at(i, j-1) + image.at(i, j+kernel.w-1);
      }
      res_1.set(i, j, sum);
    }
  }
  sum = 0;
  for (int j = 0; j < res_1.w; j++) {
    for (int i = 0; i + kernel.h <= res_1.h; i++) {
      if (i == 0) {
        sum = 0;
        for (int ii = 0; ii < i + kernel.h; ii++) {
          sum += res_1.at(ii, j);
        }
      }
      else {
        sum = sum - res_1.at(i-1, j) + res_1.at(i+kernel.h-1, j);
      }
      res_2.set(i, j, sum/(kernel.h*kernel.w));
    }
  }
  return res_2;
}

// trivial, but parallel
matrix trivial_parallel(matrix& image, matrix& kernel) {
  assert(kernel <= image);
  matrix res(image.h-kernel.h+1, image.w-kernel.w+1, true);
  #pragma omp parallel for
  for (int i = 0; i <= image.h - kernel.h; i++) {
    for (int j = 0; j + kernel.w <= image.w; j++) {
      double avg = 0;
      for (int ii = i; ii < i + kernel.h; ii++) {
        for (int jj = j; jj < j + kernel.w; jj++) {
          avg += image.at(ii, jj);
        }
      }
      avg /= (kernel.h*kernel.w);
      res.set(i, j, avg);
    }
  }
  return res;
}

// memorize, and parallel
matrix memorize_parallel(matrix& image, matrix& kernel) {
  assert(kernel <= image);
  matrix res_1(image.h, image.w-kernel.w+1, true);
  matrix res_2(image.h-kernel.h+1, image.w-kernel.w+1, true);
  #pragma omp parallel for
  for (int i = 0; i < image.h; i++) {
    double sum;
    for (int j = 0; j + kernel.w <= image.w; j++) {
      if (j == 0) {
        sum = 0;
        for (int jj = 0; jj < j + kernel.w; jj++) {
          sum += image.at(i, jj);
        }
      }
      else {
        sum = sum - image.at(i, j-1) + image.at(i, j+kernel.w-1);
      }
      res_1.set(i, j, sum);
    }
  }
  #pragma omp parallel for
  for (int j = 0; j < res_1.w; j++) {
    double sum;
    for (int i = 0; i + kernel.h <= res_1.h; i++) {
      if (i == 0) {
        sum = 0;
        for (int ii = 0; ii < i + kernel.h; ii++) {
          sum += res_1.at(ii, j);
        }
      }
      else {
        sum = sum - res_1.at(i-1, j) + res_1.at(i+kernel.h-1, j);
      }
      res_2.set(i, j, sum/(kernel.h*kernel.w));
    }
  }
  return res_2;
}



int main() {

  matrix image(1000, 1000, true);
  image.random_init();
  matrix kernel(10, 10);

  auto startTime = std::chrono::high_resolution_clock::now();
  matrix res1 = trivial(image, kernel);
  auto endTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by trivial: " << duration << " microseconds." << std::endl;

  startTime = std::chrono::high_resolution_clock::now();
  matrix res2 = memoize(image, kernel);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by memoize: " << duration << " microseconds." << std::endl;

  startTime = std::chrono::high_resolution_clock::now();
  matrix res3 = trivial_parallel(image, kernel);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by trivial_parallel: " << duration << " microseconds." << std::endl;

  startTime = std::chrono::high_resolution_clock::now();
  matrix res4 = memorize_parallel(image, kernel);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  std::cout << "Time taken by memorize_parallel: " << duration << " microseconds." << std::endl;


  // check
  assert(res1 == res2);
  assert(res1 == res3);
  assert(res1 == res4);
}