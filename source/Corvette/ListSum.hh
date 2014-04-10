#ifndef __LISTSUM__
#define __LISTSUM__

#include <vector>


template<typename Float, typename DoubleFloat>
class ListSum {
  typedef std::vector<Float> Array;
  typedef typename Array::size_type ArrayIndex;
  DoubleFloat Sum;
  std::vector<Float> Data;

  enum {ResetPeriod = 1000000};   // Number of updates before a rebuild. It is defined a constant here with the enum trick.
  long long ResetCountdown;       // Count down before next resetting of energies. This is only used to fix accumulating floating point errors for the energies

  inline void reset() {
    ResetCountdown = ResetPeriod;
    Sum = 0;
    for (ArrayIndex i = 0; i < Data.size(); ++i)
      Sum += Data[i];
  }

public:
  ListSum(const ArrayIndex& nterms) : Sum(0), Data(nterms, 0), ResetCountdown(ResetPeriod) {}

  ListSum() : Sum(0), ResetCountdown(ResetPeriod) {}

  inline void resize(const ArrayIndex& nterms) {
    Data.resize(nterms, 0);
  }

  inline void update(const ArrayIndex& i, const Float& v) {
    Sum += v - Data[i];
    Data[i] = v;

    if (--ResetCountdown == 0)
      reset();
  }

  inline const DoubleFloat& value() const {
    return Sum;
  }

};

#endif // __LISTSUMCONTAINER__