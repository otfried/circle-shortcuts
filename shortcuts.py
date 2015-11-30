#
# Calculations for "Shortcuts for the circle"
# 
# run with Python 2
#

from __future__ import division, print_function
from math import *

# Lemma numbers

lem_antipodal = 9
lem_shortlong = 10

# alpha(a)
def shortcutAngle(a):
  return 2 * asin(a/2)

def aFromAlpha(alpha):
  return 2 * sin(alpha/2)
  
def delta(a):
  return asin(a/2) - a/2

# compute a s.t. delta(a) = d
# by binary search  
def aFromDelta(d):
  a0 = 0
  a1 = 2.0
  while a1 - a0 > 0.00000001:
    am = 0.5 * (a0 + a1)
    dm = delta(am)
    if dm > d:
      a1 = am
    else:
      a0 = am
  return a0

# compute a* by binary search
def aStarKFromK(k):
  a0 = 0
  a1 = 2.0
  target = (k-1) * pi / k
  while a1 - a0 > 0.00000001:
    am = 0.5 * (a0 + a1)
    tm = am + delta(am)
    if tm > target:
      a1 = am
    else:
      a0 = am
  return a0

def deepFromA(a):
  alpha = shortcutAngle(a)
  x0 = 0
  x1 = aFromAlpha(alpha/2)
  while x1 - x0 > 0.00000001:
    xm = 0.5 * (x0 + x1)
    arc = shortcutAngle(xm)
    if xm + a > alpha - arc: 
      x1 = xm
    else:
      x0 = xm
  return alpha - 2*arc

def deepumbra():
  print("Section 2:")
  print("==========\n")
  print("Computing the length of the deep umbra:")
  for a in [0.1, 0.2, 0.5, 1.0, 1.5, 1.8, 1.9, 2.0]:
    print("a = %g,  x = %1.8f" % (a, deepFromA(a)))

def uptofive1():
  print()
  print("Section 3.1: No combinations")
  print("============================\n")
  print("Here we compute Table 1, with a*, d*, and pi - d*")
  print("We already include the values for k = 6")
  for k in range(2, 6):
    ask = aStarKFromK(k)
    print("%d  %1.4f %1.4f %1.4f" % (k, ask, delta(ask), pi - delta(ask)))
  print("%d  %1.4f %1.4f %1.4f" % (6, 2.0, delta(2.0), pi - delta(2.0)))

def uptofive2():
  print("\nSection 3.2: Antipodal pair combinations")
  print(  "========================================\n")
  print("Computing mu_k")
  for k in range(2, 6):
    ask = aStarKFromK(k)
    dsk = delta(ask)
    muk = aFromDelta(dsk/2)
    print("%d  %1.4f" % (k, muk))
  print("%d  %1.4f" % (6, aFromDelta(delta(2)/2)))
  antipodal_k3()
  lemma_shortlong()
  antipodal_k456()

def antipodal_k3():
  print("\nLemma %d for k = 3\n" % lem_antipodal)
  a = aStarKFromK(3)
  ds = delta(a)
  mu = aFromDelta(ds/2)
  print("a* = %1.4f, d* = %1.4f, mu = %1.4f" % (a, ds, mu))
  print("(i) (pi - d*)/2 = %1.4f" % ((pi - ds)/2))
  print("(v) delta(1.45) = %1.4f" % delta(1.45))
  print("(vi) delta(pi/2) = %1.4f" % delta(pi/2))
  print("(vii) a such that delta(a) = 0.06 = %1.4f" % aFromDelta(0.06))
  print("0.54 is %g degrees" % degrees(0.54))

# This computation leads to the same results
def shortLong_alternative(ds):
  # We need to solve delta(x) + delta(y) = ds/2,
  # where x^2 + y^2 = 4
  x0 = 0
  x1 = sqrt(2)
  target = ds/2
  while x1 - x0 > 0.00000001:
    x = 0.5 * (x0 + x1)
    tm = delta(x) + delta(sqrt(4 - x*x))
    if tm > target:
      x0 = x
    else:
      x1 = x
  return x0, sqrt(4 - x0*x0)

def shortLong(ds):
  # We need to solve delta(x) + delta(y) = ds/2,
  # while x + y = pi - ds
  # so y = pi - ds - x
  x0 = 0
  x1 = (pi - ds) / 2
  target = ds/2
  while x1 - x0 > 0.00000001:
    x = 0.5 * (x0 + x1)
    y = min(pi - ds - x, 2)
    tm = delta(x) + delta(y)
    if tm > target:
      x0 = x
    else:
      x1 = x
  return x0, pi - ds - x0

def lemma_shortlong():
  print("\nLemma %d: sigma_k and lambda_k\n" % lem_shortlong)
  print("Here we compute Table 2")
  print("k & a* & d* & mu_k & sigma_k & lambda_k")
  for k in range(4, 7):
    a = aStarKFromK(k)
    ds = delta(a)
    mu = aFromDelta(ds/2)
    sk, lk = shortLong(ds)
    print("%d & %1.4f & %1.4f & %1.4f & %1.4f & %1.4f \\\\" %
          (k, a, ds, mu, sk, lk))
  print("\nCalculations in the lemma:\n")
  for k in range(4, 7):
    a = aStarKFromK(k)
    ds = delta(a)
    sk, lk = shortLong(ds)
    print(" k=%d: 2 * delta(lambda_k) = %1.4f, d*/2 = %1.4f" %
          (k, 2*delta(lk), ds/2))

def antipodal_k456():
  print("\nLemma %d for k in {4, 5, 6}\n" % lem_antipodal)
  print("Showing that l = 1:\n")
  for k in range(4, 7):
    a = aStarKFromK(k)
    ds = delta(a)
    mu = aFromDelta(ds/2)
    sk, lk = shortLong(ds)
    dh = ds - 2 * (k - 3) * delta(sk)
    w = pi - lk - dh
    check1 = (k-2) * w
    check2 = (k-1) * w
    form = "%d %1.4f %1.4f %1.4f"
    print("k = %d, d* = %1.4f, sigma = %1.4f, lambda = %1.4f" %
          (k, ds, sk, lk))
    print("  delta^ = %1.4f,  w = (pi - lambda - d^) = %1.4f" % (dh, w))
    print("  (k-2) * w = %1.4f < pi" % ((k-2) * w))
  print("\nThe final contradiction of Lemma 6:\n")
  for k in range(4, 7):
    a = aStarKFromK(k)
    ds = delta(a)
    mu = aFromDelta(ds/2)
    sk, lk = shortLong(ds)
    dh = ds - 2 * (k - 3) * delta(sk)
    ah = aFromDelta(dh)
    w = pi - ah - dh
    print("k = %d, delta^ = %1.4f, a^ = %1.4f, w = pi - a^ - delta^ = %1.4f" %
          (k, dh, ah, w))
    print("   (k-1) * w = %1.4f < pi" % ((k-1) * w))

def uptofive3():
  print("\nSection 3.3: Combinations allowed")
  print(  "=================================\n")
  for k in range(3, 6):
    ask = aStarKFromK(k)
    dsk = delta(ask)
    muk = aFromDelta(dsk/2)
    w = pi - muk - dsk
    print("k=%d mu=%1.4f w = pi - mu - d* = %1.4f (k-1)*w = %1.4f < pi" % 
          (k, muk, w, (k-1)*w))

print("\\begin{verbatim}")
print("Source code at: http://github.com/otfried/circle-shortcuts\n")
deepumbra()
uptofive1()
uptofive2()
uptofive3()
print("\\end{verbatim}")

