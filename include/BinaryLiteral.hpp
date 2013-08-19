#ifndef __BINARY_LITERAL_HPP__
#define __BINARY_LITERAL_HPP__

//===============================================================

//
// binary_literal_impl represents a compile-time function that
// computes the unsigned long long int from a list of characters
// Digits that MUST be composed of '0' or '1'...
//
template <char... Digits>
struct binary_literal_impl;

// If the next digit is zero, then compute the rest...
template <char... Digits>
struct binary_literal_impl<'0', Digits...>
{
  static constexpr unsigned long long to_ulonglong()
  {
    return binary_literal_impl<Digits...>::to_ulonglong();
  }
};

// If the next digit is one, then shift 1 and compute the rest...
template <char... Digits>
struct binary_literal_impl<'1', Digits...>
{
  static constexpr unsigned long long to_ulonglong()
  {
    return (1UL << sizeof...(Digits))
      | binary_literal_impl<Digits...>::to_ulonglong();
  }
};

// Base case: No digits, so return 0...
template <>
struct binary_literal_impl<>
{
  static constexpr unsigned long long to_ulonglong()
  {
    return 0;
  }
};

//===============================================================

template <char... Digits>
constexpr unsigned long long operator "" _binary()
{
  return binary_literal_impl<Digits...>::to_ulonglong();
}

//===============================================================

#endif // __BINARY_LITERAL_HPP__