


# base/intfuncs.jl provides these number theoretic functions:
# * gcd(x,y)
# * isprime(x)
# * powermod(x, a, m) - integer power under the modulus m
# * invmod(x, m) - returns the reciprocal of an integer mod m


# Primes.jl provides:
# * primes
# * nextprime(x+1)
# * prevprime(x-1)
# * totient(x)
# * factor(Vector, x)
using Primes


function ntt(xvec, root, mod)
  if length(xvec) >= mod
    Base.error("mod must be greater than length of input vector");
  end

  if !(1 <= root < mod)
    Base.error("!(1 <= root < mod)");
  end

  # yvec = Vector{UInt64}(undef, length(xvec));
  yvec = Vector{Complex{UInt64}}(undef,length(xvec));

  for ii in 1:length(xvec)
    # y = zero(UInt64);
    y = Complex{UInt64}(0 + 0im);
    for (jj, x) in enumerate(xvec)
      println("loop($ii, $jj)");
      # y = (y + x * powermod(root, UInt64(ii*jj), mod)) % mod;
      y = (y + x * powermod(root, UInt64(ii*jj), mod));
      y = (y.re % mod) + (y.im % mod)im;
    end
    yvec[ii] = y;
  end

  return yvec;
end



"""
Generates a table of prime numbers such that:
 N is a prime number
 2 is a primitive (N-1)th root of unity.
Also computes NTT vector lengths n for which the prime number N is valid.
Choose an arbitrary integer k
 k >= 1
 N = k*n + 1
All n that satisfy this equation are valid vector lengths.

Also computes the maximum input value m for a length n vector convolution
that would avoid overflow in the worst case.
 N >= m^2 * n + 1

@param [in] N Prime number.

Example:
Prime: 2147483867
 2 is 2147483866 primitive root of unity.
 ... ignore a bunch and get to the good stuff ...
 25786.0  max-value: 288.0
 15142.0  max-value: 376.0
 12893.0  max-value: 408.0
 7571.0  max-value: 532.0
 2486.0  max-value: 929.0
 1474.0  max-value: 1207.0
 1243.0  max-value: 1314.0
 737.0  max-value: 1706.0
 226.0  max-value: 3082.0
 134.0  max-value: 4003.0
 113.0  max-value: 4359.0
 67.0  max-value: 5661.0
 22.0  max-value: 9879.0

 Prime: 2147484541
 2 is 2147484540 primitive root of unity.
 1020.0  max-value: 1450.0
 510.0  max-value: 2052.0
 340.0  max-value: 2513.0
 255.0  max-value: 2901.0
 204.0  max-value: 3244.0
 170.0  max-value: 3554.0
 102.0  max-value: 4588.0
 85.0  max-value: 5026.0
 68.0  max-value: 5619.0
 60.0  max-value: 5982.0
 51.0  max-value: 6489.0
 34.0  max-value: 7947.0
 30.0  max-value: 8460.0
 20.0  max-value: 10362.0
 17.0  max-value: 11239.0

 Prime: 2147484949
  2 is 2147484948 primitive root of unity.
  2919.0  max-value: 857.0
  2604.0  max-value: 908.0
  2443.0  max-value: 937.0
  2363.0  max-value: 953.0
  2108.0  max-value: 1009.0
  2094.0  max-value: 1012.0
  1946.0  max-value: 1050.0
  1668.0  max-value: 1134.0
  1581.0  max-value: 1165.0
  1428.0  max-value: 1226.0
  1396.0  max-value: 1240.0
  1302.0  max-value: 1284.0
  1054.0  max-value: 1427.0
  1047.0  max-value: 1432.0
  973.0  max-value: 1485.0
  868.0  max-value: 1572.0
  834.0  max-value: 1604.0
  714.0  max-value: 1734.0
  698.0  max-value: 1754.0
  651.0  max-value: 1816.0
  556.0  max-value: 1965.0
  527.0  max-value: 2018.0
  476.0  max-value: 2124.0
  434.0  max-value: 2224.0
  417.0  max-value: 2269.0
  372.0  max-value: 2402.0
  357.0  max-value: 2452.0
  349.0  max-value: 2480.0
  278.0  max-value: 2779.0
  238.0  max-value: 3003.0
  217.0  max-value: 3145.0
  204.0  max-value: 3244.0
  186.0  max-value: 3397.0
  139.0  max-value: 3930.0
  124.0  max-value: 4161.0
  119.0  max-value: 4248.0
  102.0  max-value: 4588.0


"""
function ntt_prime_describe(N)
  println("Prime: $N");
  N1 = N-1;
  if totient(N) != N1
    println(" totient($N) != $N-1");
    return nothing;
  end

  if !is_generator(2, N1, N)
    println(" 2 is not a generator for this prime");
    return nothing;
  end

  println(" 2 is $N1 primitive root of unity.");

  k = 1;
  n = N1;
  while n > 15 # some small length vector that I will never do an NTT on
    if n < 3000 && isinteger(n)
      in_max_value = floor(sqrt(N1/n));
      println(" $n  max-value: $in_max_value")
    end
    k = k + 1;
    n = N1/k;
  end
end


"""
Returns an arbitrary generator of the multiplicative integer group modulo m.
An answer must exist if m is prime.

@param [in] totient Must equal the Euler phi function.
@param [in] m Modulo.
"""
function find_generator(totient, m)
  if !(1 <= totient < m)
    Base.error("!(1 <= totient < m)");
  end

  println("find_generator($totient,$m): begin");
  for ii in 1:(m-1)
    if is_generator(ii, totient, m)
      println("  find_generator($totient,$m): $ii");
      return ii;
    end
  end
  println("find_generator($totient,$m): done");

  Base.error("No generator exists.")
end


"""

@param [in] n Degree.
@param [in] totient Must be a multiple of degree.
@param [in] m Modulo.
"""
function find_primitive_root(n, totient, m)
  if !(1 <= n <= totient < m)
    Base.error("!(1 <= n <= totient < m)");
  end

  if totient % n != 0
    Base.error("totient % n != 0");
  end

  g = find_generator(totient, m);
  root = powermod(g, Int64(floor(totient/n)), m);

  if !(0 <= root < m)
    Base.error("!(0 <= root < m)");
  end

  return root;
end



"""
Tests whether x generates the multiplicative group of integers modulo m.
If x is a generator then it produces all integers in the range [1 mod).
"""
function is_generator(x, totient, m)
  if !(0 <= x < m)
    Base.error("!(0 <= x < mod)");
  end

  if !(1 <= totient < m)
    Base.error("!(1 <= totient < m)");
  end

  pfs = unique_prime_factors(totient);
  return (powermod(x, totient, m) == 1) && all([powermod(x, Int64(floor(totient/pf)), m) != 1 for pf in pfs]);
end


"""
Tests whether x is a primitive n-th degree root of unity modulo m.
  x^n % m = 1
  x^k % m != 1 for 1 <= k < n

@param [in] x Value to be tested.
@param [in] n Degree.
@param [in] m Modulus.
"""
function is_primitive_root(x, n, m)::Bool
  if !(0 <= x < m)
    Base.error("!(0 <= x < m)");
  end

  if !(1 <= n < m)
    Base.error("!(1 <= n < m)");
  end

  pfs = unique_prime_factors(n);
  return (powermod(x, n, m) == 1) && all([powermod(x, Int64(floor(n/pf)), m) != 1 for pf in pfs]);
end


# # Use the julia builtin invmod() function.
# """
# Calculates the multiplicative inverse of x modulo m.
# The modular multiplicative inverse is an integer y such that:
#   0 <= y < m
#   (y*x) % m = 1
# The multiplicative inverse exists if and only if x and m are relatively prime:
#   gcd(x, m) = 1.
#
# This function assumes m is prime and uses "Fermat's Little theorem" to solve.
# If x and m are relatively prime, then modulo inverse is:
#   x^(m-2) % m
#
# """
# function invmod(x, m)
#   if gcd(x, m) != 1
#     Base.error("gcd($x,$m) != 1");
#   end
#
#   return powermod(x, m-2, m);
# end  # invmod


"""
Returns a list of unique prime factors in ascending order.
Example: Prime factors of 60 are:
  2^2 ⋅ 3 ⋅ 5
The unique prime factors are:
  2, 3, 5

@param [in] x Numeric value. Will be rounded to nearest integer.
"""
function unique_prime_factors(x::T)::Vector{T} where {T <: Integer}
  return unique(factor(Vector, x));
end

function unique_prime_factors(x::T)::Vector{Int64} where {T <: Number}
  return unique_prime_factors(Int64(round(x)));
end
