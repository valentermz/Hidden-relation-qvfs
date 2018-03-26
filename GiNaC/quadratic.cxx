/** Find the “hidden” relation described in https://arxiv.org/abs/1705.06340
 *
 * Copyright (C) 2017, 2018 Yury G. Kudryashov
 *
 * This file is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 2 of the License or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this
 * program. If not, see <http://www.gnu.org/licenses/>
 *
 * Note that GiNaC is GPL v2 only, so if you want to use GPL v3 or later, you need
 * to port the code to some other library.
 */
#include <ginac/ginac.h>
#include <array>
#include <cassert>
#include <utility>

using namespace GiNaC;
using std::cout;
using std::cerr;
using std::endl;
using std::array;
using std::pair;
using std::make_pair;

const array<pair<int, int>, 3> p {make_pair(0, 0), make_pair(1, 0), make_pair(0, 1)};

/// A quadratic polynomial vanishing at @c p
class QuadPoly
{
private:
  ex m_x, m_y, m_xy;
public:
  QuadPoly(ex x, ex y, ex xy):
    m_x {x}, m_y {y}, m_xy {xy}
  { }

  ex operator()(const ex &x, const ex &y) const
  {
    return m_x * x + m_y * y - m_x * x * x + m_xy * x * y - m_y * y * y;
  }

  ex at_infty(const ex &w) const
  {
    symbol z{"z"};
    return ((*this)(1 / z, w / z) * z * z).normal().subs({{z, 0}});
  }
};

void test_quadpoly()
{
  symbol cx, cy, cxy, x, y;
  QuadPoly P{cx, cy, cxy};
  for(const auto &z: p)
    assert(P(z.first, z.second) == 0);
}

/// A quadratic vector field vanishing at @c p
class QuadVF
{
private:
  QuadPoly m_P, m_Q;
  exmap m_sub;
public:
  QuadVF(QuadPoly P, QuadPoly Q):
    m_P(P), m_Q(Q)
  { }

  void add_sub(const exmap &other)
  {
    m_sub.insert(other.begin(), other.end());
  }

  const exmap &sub() const
  { return m_sub; }

  matrix J(size_t i) const
  {
    symbol x, y;
    ex P = m_P(x, y);
    ex Q = m_Q(x, y);
    exmap pt {{x, p[i].first}, {y, p[i].second}};
    matrix res
      {
	{P.diff(x).subs(pt), P.diff(y).subs(pt)},
	  {Q.diff(x).subs(pt), Q.diff(y).subs(pt)}
      };
    return ex_to<matrix>(res.subs(m_sub));
  }

  ex det(size_t i) const
  {
    return J(i).determinant().subs(m_sub).normal();
  }

  ex prod_l_infty() const
  {
    symbol w {"w"};
    ex F = m_P.at_infty(w);
    ex G = w * F - m_Q.at_infty(w);
    ex Gp = G.diff(w);
    ex res = (resultant(F, G, w) / resultant(Gp, G, w)).subs(m_sub).numer_denom();
    return (res.op(0).expand().subs(m_sub) / res.op(1).expand().subs(m_sub));
  }

  static QuadVF from_spectra(const array<ex, 3> &t, const array<ex, 3> &d, const symbol &sqD);
};

ex lsolve1(const ex &eq, const symbol &s)
{
  return -eq.coeff(s, 0) / eq.coeff(s, 1);
}

QuadVF QuadVF::from_spectra(const array<ex, 3> &t, const array<ex, 3> &d, const symbol &sqD)
{
  // Introduce some symbols for temporary unknowns.
  symbol px {"px"}, py {"py"}, qx {"qx"};
  // Take polynomials P and Q with prescribed traces
  ex qy = t[0] - px;
  QuadPoly P {px, py, t[2] + qy - px};
  QuadPoly Q {qx, qy, t[1] + px - qy};
  // Assert that our formulas are correct
  QuadVF v {P, Q};
  for(size_t i = 0; i < 3; ++ i)
    assert(v.J(i).trace() == t[i]);
  // Now, take determinants into account
  array<ex, 3> eqs_det;
  for(size_t i = 0; i < 3; ++ i)
      eqs_det[i] = v.J(i).determinant() - d[i];
  ex v_py = lsolve1(eqs_det[0] + eqs_det[2], py);
  ex v_qx = lsolve1(eqs_det[0], qx).subs({{py, v_py}});
  ex eq = eqs_det[1].subs({{qx, v_qx}, {py, v_py}}).normal().numer();
  ex a = eq.coeff(px, 2);
  ex b = eq.coeff(px, 1);
  ex c = eq.coeff(px, 0);
  ex D = b * b - 4 * a * c;
  exmap s;
  for(size_t i = 1; i < 4; ++ i) {
    s[pow(sqD, 2 * i)] = pow(D, i);
    s[pow(sqD, 2 * i + 1)] = sqD * pow(D, i);
  }
  s[px] = ((-b + sqD) / (2 * a)).normal();
  s[py] = v_py.subs(s).normal().subs(s).normal();
  s[qx] = v_qx.subs(s).normal().subs(s).normal();
  v.add_sub(s);
  for(const auto &i: s)
    cerr << ex(i.first == i.second) << endl;
  for(size_t i = 0; i < 3; ++ i)
    assert(v.det(i) == d[i]);
  return v;
}

/** Find A such that A[0] * denom^2 + A[1] * numer * denom + A[2] * numer^2 = 0.
 *
 * Here numer and denom are linear in x, and x^2=D.
 */
std::array<ex, 3> find_equation(const ex &numer, const ex &denom, const ex &D, const symbol &x)
{
  assert(numer.degree(x) == 1);
  assert(denom.degree(x) == 1);
  ex An = numer.coeff(x, 1), Bn = numer.coeff(x, 0), Ad = denom.coeff(x, 1), Bd = denom.coeff(x, 0);
  return {
    Bd * Bd - Ad * Ad * D,
      -2 * (Bn * Bd - An * Ad * D),
      Bn * Bn - An * An * D
      };
}

void test_find_equation()
{
  std::array<symbol, 4> c {symbol {"a"}, symbol {"b"}, symbol {"c"}, symbol {"d"}};
  symbol D {"D"}, x {"√D"};
  ex n = c[0] * x + c[1];
  ex d = c[2] * x + c[3];
  auto res = find_equation(n, d, D, x);
  exmap s;
  s[x * x] = D;
  ex eq = (res[0] * n * n + res[1] * n * d + res[2] * d * d).expand().subs(s).expand();
  cerr << eq << endl;
  assert(eq == 0);
}


int main(int argc, char *argv[])
{
  test_find_equation();
  array<ex, 3> t {symbol{"t0"}, symbol{"t1"}, symbol{"t2"}};
  array<ex, 3> d {symbol{"d0"}, symbol{"d1"}, symbol{"d2"}};
  symbol sqD {"√D"};
  QuadVF v = QuadVF::from_spectra(t, d, sqD);
  ex D = v.sub().at(sqD * sqD);
  cerr << "Constructed a vector field" << endl;
  ex Lambda_nd = v.prod_l_infty().numer_denom();
  cerr << Lambda_nd << endl;
  auto H = find_equation(Lambda_nd.op(0), Lambda_nd.op(1), D, sqD);
  symbol Lambda("Λ");
  ex f = gcd(H[0], gcd(H[1], H[2]));
  cout << "H:" << endl;
  for(auto &h: H) {
    h = factor((h / f).normal());
    cout << "\t" << h << endl;
  }

  return 0;
}
