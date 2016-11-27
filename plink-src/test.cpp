#include <iostream>

using namespace std;


class A {
public:
  virtual void fitModel() = 0;
  void generalFunction() { cout << "General function (A)\n"; }
};

class B : public A {
public:
  void fitModel() { cout << "Fit model (B)\n"; }
};

class C : public A {
public:
  void fitModel() { cout << "Fit model (C)\n"; }
};

int main() {

  A * a;

  if (false)
    {
      B * b = new B;
      a = b;
    }
  else
    {
      C * c = new C;
      a = c;
    }	
  
  a->generalFunction();
  a->fitModel();
  
  exit(0);
}
