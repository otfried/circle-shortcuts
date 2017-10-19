#
# Calculations for "Shortcuts for the circle"
# 
# run with Python 2 or 3
#

from __future__ import division, print_function
from math import *

# Lemma numbers

lem_deep = 4
lem_antipodal = 10
lem_upto6combined = 12
thm_6to7 = 14
lem_Sp6 = 16
lem_s2dht = 17
lem_s2longvi = 18
lem_thetaset = 19
lem_bound_combining_edges = 20

# precision for binary search
precision = 1e-12

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
  while a1 - a0 > precision:
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
  while a1 - a0 > precision:
    am = 0.5 * (a0 + a1)
    tm = am + delta(am)
    if tm > target:
      a1 = am
    else:
      a0 = am
  return a0

def umbra():
  print("Lemma %d:" % lem_deep)
  print("=========\n")
  print("delta(delta(2)) = %1.5f" % delta(delta(2)))

def uptofive1():
  print("\nSection 3.1:")
  print(  "===========\n")
  print("Table 1:")
  print("k  a*     d*     pi - d*")
  for k in range(2, 6):
    ask = aStarKFromK(k)
    print("%d  %1.4f %1.4f %1.4f" % (k, ask, delta(ask), pi - delta(ask)))
  print("%d  %1.4f %1.4f %1.4f" % (6, 2.0, delta(2.0), pi - delta(2.0)))
  print("\nComputing mu_k:")
  for k in range(2, 6):
    ask = aStarKFromK(k)
    dsk = delta(ask)
    muk = aFromDelta(dsk/2)
    print("mu_%d = %1.4f" % (k, muk))
  print("mu_6 = %1.4f" % (aFromDelta(delta(2)/2)))


def uptofive2():
  print("\nLemma %d:" % lem_antipodal)
  print(  "=========\n")
  antipodal_k3()
  table2()
  antipodal_k456()

def antipodal_k3():
  print("Lemma %d for k = 3\n" % lem_antipodal)
  a = aStarKFromK(3)
  ds = delta(a)
  mu = aFromDelta(ds/2)
  print("a* = %1.4f, d* = %1.4f, mu = %1.4f" % (a, ds, mu))
  print("(i) (pi - d*)/2 = %1.4f" % ((pi - ds)/2))
  print("(v) delta(1.45) = %1.4f" % delta(1.45))
  print("(vi) delta(pi/2) = %1.4f" % delta(pi/2))
  print("(vii) a such that delta(a) = 0.06 = %1.4f" % aFromDelta(0.06))

def shortLong(ds):
  # We need to solve delta(x) + delta(pi - ds - x) = ds/2,
  # for x in [pi - ds - 2, (pi-ds)/2]
  x0 = pi - ds - 2
  x1 = (pi - ds) / 2
  target = ds/2
  while x1 - x0 > precision:
    x = 0.5 * (x0 + x1)
    tm = delta(x) + delta(pi - ds - x)
    if tm > target:
      x0 = x
    else:
      x1 = x
  return x0, pi - ds - x0

def table2():
  print("\nTable 2:")
  print("k & a*     & d*     & mu_k   & sigma_k & lambda_k")
  for k in range(4, 7):
    a = aStarKFromK(k)
    ds = delta(a)
    mu = aFromDelta(ds/2)
    sk, lk = shortLong(ds)
    print("%d & %1.4f & %1.4f & %1.4f & %1.4f  & %1.4f \\\\" %
          (k, a, ds, mu, sk, lk))

def antipodal_k456():
  print("\nLemma %d for k in {4, 5, 6}:\n" % lem_antipodal)
  print("Showing that l = 1:")
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
    print("k = %d:" % k)
    print("  d* = %1.4f, sigma = %1.4f, lambda = %1.4f" % (ds, sk, lk))
    print("  delta^ = %1.4f,  w = (pi - lambda - d^) = %1.4f" % (dh, w))
    print("  (k-2) w = %1.4f < pi" % ((k-2) * w))
  print("\nThe final contradiction of Lemma %d:" % lem_antipodal)
  for k in range(4, 7):
    a = aStarKFromK(k)
    ds = delta(a)
    mu = aFromDelta(ds/2)
    sk, lk = shortLong(ds)
    dh = ds - 2 * (k - 3) * delta(sk)
    ah = aFromDelta(dh)
    w = pi - ah - dh
    print("k = %d:" % k)
    print("  delta^ = %1.4f, a^ = %1.4f, w = pi - a^ - delta^ = %1.4f" %
          (dh, ah, w))
    print("  (k-1) w = %1.4f < pi" % ((k-1) * w))

def uptofive3():
  print("\nProof of Lemma %d:" % lem_upto6combined)
  print(  "==================\n")
  for k in range(3, 6):
    ask = aStarKFromK(k)
    dsk = delta(ask)
    muk = aFromDelta(dsk/2)
    w = pi - muk - dsk
    print("k=%d: mu=%1.4f, w = pi - mu - d* = %1.4f => (k-1)w = %1.4f < pi" % 
          (k, muk, w, (k-1)*w))

def six():
  print("\nProof of Theorem %d for k = 6:" % thm_6to7)
  print(  "==============================\n")
  ask = aStarKFromK(6)
  dsk = delta(ask)
  muk = aFromDelta(dsk/2)
  w = pi - muk - dsk
  print("mu=%1.4f, d*=%1.4f, w = pi - mu - d* = %1.4f" % (muk, dsk, w))
  print(" 4 w = %1.4f < pi" % (4 * w))
  print(" pi - d* + 5 w = %1.4f < 2pi" % (pi - dsk + 5 * w))
  print(" delta(pi - d* - mu) = %1.4f" % (delta(w)))
  dht = dsk - 0.016
  aht = aFromDelta(dht)
  print(" d^ = %1.4f, a^ = %1.4f" % (dht, aht))
  print(" 5 (pi - a^ - d^) = %1.4f < pi" % (5 * (pi - aht - dht)))

def seven():
  print("\nSeven shortcuts, Proof of Lemma %d:" % lem_Sp6)
  print(  "===================================\n")
  ds = pi/2 - 1
  s6, l6 = shortLong(ds)
  dht4 = ds - 8 * delta(s6)
  print("d* = %1.4f, sigma6 = %1.4f, delta(sigma6) = %1.4f" %
        (ds, s6, delta(s6)))
  print("d^(4) = %1.4f, 2 delta(1.849) = %1.4f" % (dht4, 2 * delta(1.849)))
  w = pi - 1.849 - dht4
  print("  w = pi - 1.849 - d^(4) = %1.4f" % w)
  print("  4 * w = %1.4f < pi" % (4 * w))
  dht2 = ds - 4 * delta(s6)
  ds5 = delta(aStarKFromK(5))
  print("d^(2) = %1.4f > %1.4f = d5*" % (dht2, ds5))
  dht = ds - 2 * delta(s6)
  print("d^(1) = %1.4f" % dht)
  w = pi - l6 - dht
  print("  w = pi - lambda6 - d^(1) = %1.4f" % w)
  print("  4 * w = %1.4f < pi" % (4 * w))
  #
  print("\nSeven shortcuts, Proof of Lemma %d:" % lem_s2dht)
  print(  "===================================\n")
  v = 1.7
  print("1.7 + lambda6 = %1.4f > pi" % (v+l6))
  print("delta(1.999) = %1.4f < 0.54" % delta(1.999))
  w1 = pi - 1.999 - 0.54
  w2 = pi - l6 - 0.54
  w3 = pi - v - 0.54
  print("  w1 = pi - 1.999 - 0.54 = %1.4f" % w1)
  print("  w2 = pi - lambda6 - 0.54 = %1.4f" % w2)
  print("  w3 = pi - 1.7 - 0.54 = %1.4f" % (w3))
  print("  w3 + 2 * w2 + 6 * w1 = %1.4f < 2pi" % (w3+2*w2+6*w1))
  print("  2 * w2 + 8 * w1 = %1.4f < 2pi" % (2 * w2 + 8 * w1))
  # 
  print("\nSeven shortcuts, Proof of Lemma %d:" % lem_s2longvi)
  print(  "===================================\n")
  dht = ds - 2 * delta(s6)
  print("d^ = d* - 2 * delta(sigma6) = %1.4f" % dht)
  s3 = 0.8 * pi - dht
  print("s3 <= 0.8 * pi - d^ = %1.4f" % s3)
  A1 = 4 * delta(l6) * (pi - l6 - dht)
  A2 = 4 * delta(s3) * (pi - s3 - dht)
  A3 = 4 * dht * (pi - 1.999 - dht)
  print("A(lambda6, d^) = 4 * delta(lambda6) * (pi - lambda6 - d^) = %1.4f"
        % A1)
  print("A(s3, d^) <= 4 * delta(%1.4f) * (pi - %1.4f - d^) = %1.4f" %
        (s3,s3,A2))
  print("delta(1.999) = %1.5f < d^" % delta(1.999))
  print("A(1.999, d^) < 4 * d^ * (pi - 1.999 - d^) = %1.4f" % A3)
  print("%1.4f + %1.4f + 4 * %1.4f = %1.4f < %1.4f = 4 * d^ * pi" %
        (A1, A2, A3, A1 + A2 + 4 * A3, 4 * dht * pi))
  #
  print("\nSeven shortcuts, Proof of Lemma %d:" % lem_thetaset)
  print(  "===================================\n")
  print("0.4 - d*/2 = %1.4f, 0.4 + d*/2 = %1.4f" %
        (0.4 - ds/2, 0.4 + ds/2))
  # 
  print("\nSeven shortcuts, a short shortcut exists")
  print(  "========================================\n")
  w = pi - 1.999 - ds
  print("pi - 1.999 - d* = %1.4f" % w)
  print("pi - 2 * 0.4 = %1.4f > %1.4f = 4 * %1.4f" % (pi - 2 * 0.4, 4 * w, w))
  print("0.4 + d*/2 = %1.4f > %1.4f" % (0.4 + ds/2, w))
  # 
  print("\nSeven shortcuts, Proof of Lemma %d:" % lem_bound_combining_edges)
  print(  "===================================\n")
  print("s1 + s2 <= 5 * pi - 7 * ds - 5 * lambda6 = %1.4f < 2.34" %
        (5*pi - 7*ds - 5*l6))
  print("delta(sigma6) + delta(2.34 - sigma6) = %1.4f < 0.2" %
        (delta(s6) + delta(2.34 - s6)))
  print("delta(0.83) + delta(pi/2 + 1 - 0.83) = %1.4f < 0.2" %
        (delta(0.83) + delta(pi/2 + 1 - 0.83)))
  print("delta(0.83) + delta(1.7) = %1.4f < 0.2" %
        (delta(0.83) + delta(1.7)))
  print("1.999 + sigma6 = %1.4f > %1.4f = pi - d*" %
        (1.999 + s6, pi - ds))
  # 
  print("\nSeven shortcuts, no short shortcut")
  print(  "==================================\n")
  zeta = pi/2 - 1.4
  print("zeta = pi/2 - 1.4 = %1.4f" % zeta)
  w1 = pi - l6 - ds
  print("pi - lambda6 - d* = %1.4f" % w1)
  print("2 * delta(1.949) = %1.4f < %1.4f = d* + zeta" %
        (2 * delta(1.949), ds + zeta))
  print("pi - 1.949 - d* = %1.4f" % (pi - 1.949 - ds))
  print(" 5 * 2 * 0.622 = %1.4f < 2pi" % (10 * 0.622))
  
print("\\begin{verbatim}")
print("Source code at: http://github.com/otfried/circle-shortcuts\n")
umbra()
uptofive1()
uptofive2()
uptofive3()
six()
seven()
print("\\end{verbatim}")

