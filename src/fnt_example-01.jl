




T = Int64;

N = T(4);  # Input vector length
n = T(log2(N));  # vector length power of 2

# Modulo
p = T(2^16+1);


a = T(256);

W = powermod.(a, T.((0:1:N-1) * (0:1:N-1)'), p);
Wi = powermod.(a, -T.((0:1:N-1) * (0:1:N-1)'), p);

x = T.([2, -2, 1, 0]);
# The NTT works with positive integers. Any negative integer must be made to look positive
# by adding the modulo onto it. The values can be recovered by subtracting the modulo for
# any number greater than (p-1)/2.
x = [v + Int64(v < 0) * p for v in x];


X = (W * x) .% p;
y_unscaled = (Wi * X) .% p;
y = (invmod(N, p) .* y_unscaled) .% p;



"""
Negates a number when arithmetic modulo 2^b + 1.
"""
function negate(x::Integer, b::Integer)
  mask = UInt64(2^b - 1);
  println("negate");
  println(" x: $x");
  println(" b: $b");
  println(" mask: $mask");
  return ((~x) & mask) + 2;
end

"""
When adding two b-bit integers, a b-bit integer is obtained and a possible carry bit.
The carry bit represents 2^b = -1 mod p.
To implement arithmetic modulo 2^b + 1, one must simply subtract the carry bit.
"""
function add(x::T, y::T, b::Integer) where T <: Integer
  mask = UInt64(2^b - 1);
  z = x + y;
  c = (z >> b) & one(T);
  z = z & mask;
  println("add");
  println(" x: $x");
  println(" y: $y");
  println(" c: $c");
  println(" z: $z");
  z = z - c;
  z = z & mask;
  return z;
end


"""
Subtraction is implemented as an addition by first negating the subtrahend and then
adding them.
"""
function subtract(x::T, y::T, b::Integer) where T <: Integer
  println("subtract");
  println(" x: $x");
  println(" y: $y");
  println(" b: $b");
  return add(x, negate(y, b), b);
end

"""
When the root is a power of 2, the only multiplications involved are a power of 2.
This allows multiplication to be implemented by bit shifts and addition.
@param [in] x multiplicand
@param [in] k power of 2
@param [in] b integer bit width (word size).
"""
function mul2(x::T, k::Integer, b::Integer) where T <: Integer
  println("mul2");
  mask = UInt64(2^b - 1);
  y = zero(T);
  if k >= 0
    y = x << k;
    yL = y & mask;
    yH = (y >> b) & mask;
    println(" mask: $mask");
    println(" x: $x");
    println(" k: $k");
    println(" b: $b");
    println(" y: $y");
    println(" yL: $yL");
    println(" yH: $yH");
    return subtract(yL, yH, b);
  else
    y = (x << b) >> -k;
    yL = y & mask;
    yH = (y >> b) & mask;
    println(" mask: $mask");
    println(" x: $x");
    println(" k: $k");
    println(" b: $b");
    println(" y: $y");
    println(" yL: $yL");
    println(" yH: $yH");
    return subtract(yH, yL, b);
  end

end





# END
