#include <vector>
// [[odin.dust::register]]
template <typename T, typename real_t>
real_t vaccination_schedule(size_t i, real_t daily_doses, real_t dt,
                            const T& candidates, const T& candidates_pos) {
  // Early exit in the case of no vaccination
  if (daily_doses == 0) {
    return 0;
  }
  // Early exit for our young classes:
  if (i <= 3) {
    return 0;
  }

  // 'i' comes in as base1 so we do the subtraction here so that
  // there's only one place to think about it (except for directly
  // above)
  i--;

  // Fixed priority groups for now, zero-offset indexed
  //
  // The other nice way of modelling this would be an integer matrix,
  // accepted as user-input, then loop over each row and weight.
  //
  // This way is not flexible at all, but we could switch between
  // priority groups with an enum quite efficiently.
  static const std::vector<std::vector<size_t>> priority = {
    {17, 18}, {16}, {15}, {14}, {13}, {12}, {11}, {10}, {9, 8, 7, 6, 5, 4, 3}};

  for (auto &p : priority) {
    double n = 0;
    bool exit = false;
    for (auto j : p) {
      n += candidates_pos[j];
      exit = exit || j == i;
    }
    if (exit) {
      real_t n_to_vaccinate;
      if (n < daily_doses) {
        // We won't use it all
        n_to_vaccinate = candidates_pos[i];
        // Using Inf is slightly too fast, but is in some ways better
        // than using the same value as below (with the -log) which
        // results in very soft shoulders as we hit the end of the
        // group (i.e., this branch). This is due to the way that
        // n_to_vaccinate cancels out of that equation. However, with
        // Inf we get a bit of acceleration in the final mop up stage.
        //
        // However, this will cause huge problems is we have
        // non-willing people too as they will get vaccinated!
      } else {
        // We will use it all, so share within our group
        n_to_vaccinate = daily_doses / n * candidates_pos[i];
      }
      return n_to_vaccinate * dt / candidates[i];
    } else if (n >= daily_doses) {
      // All vaccine has been used up by earlier groups
      return 0;
    } else {
      // Continue through loop
      daily_doses -= n;
    }
  }

  // Catch any unhandled early exit; this should never trigger though,
  // but the compiler will need it.
  return 0;
}
