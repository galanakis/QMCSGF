// This is a temporary function. Currently it accepts integer arguments only and returns a double.
// It needs to be generalized to complex numbers and should be moved to the Complex class file.

double Factorial(int n)
  {
    double Fact=1.0;

    while (n!=0)
      Fact*=n--;

    return Fact;
  }
